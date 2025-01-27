#!/usr/bin/env python3

import pysam
import polars as pl
import yaml

df = pl.scan_csv("${gtf}",
                    separator="\\t",
                    has_header=False,
                    new_columns=["chr", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
                    )
df = df.filter(pl.col("attributes").str.contains(r'gene_type "${gene_type}"'))
df = df.with_columns(attributes=pl.col("attributes").map_elements(lambda attributes: dict([[value.strip(r'"') for value in entry.strip().split(' ', 1)] for entry in attributes.split(';') if entry]), return_dtype=pl.Object))
df = df.with_columns(
    gene_id=pl.col("attributes").map_elements(lambda attributes: attributes.get("gene_id", None), return_dtype=pl.String),
    gene_name=pl.col("attributes").map_elements(lambda attributes: attributes.get("gene_name", None), return_dtype=pl.String),
    gene_type=pl.col("attributes").map_elements(lambda attributes: attributes.get("gene_type", None), return_dtype=pl.String),
)
df = df.drop("attributes")
df = df.filter(pl.col("gene_type") == "${gene_type}")

df = df.group_by("gene_id", "gene_name", "chr").agg(
    start=pl.col("start").min(),
    end=pl.col("end").max()
)

df = df.collect()

bam_in_path = "${bam}"
bam_in = pysam.AlignmentFile(bam_in_path, "rb")

with pysam.AlignmentFile("${prefix}.modified.bam", "wb", template=bam_in) as bam_out, pysam.AlignmentFile("${prefix}.colliders.bam", "wb", template=bam_in) as bam_colliders:
    known_reads = set()
    for row in df.select("chr", "start", "end").iter_rows():
        chr, start, end = row

        for record in bam_in.fetch(chr, start, end):
            bam_out.write(record)
            known_reads.add(record.query_name)

    for record in bam_in.fetch():
        if record.query_name in known_reads:
            bam_colliders.write(record)
        else:
            bam_out.write(record)

versions = {
    "versions": {
        "polars": pl.__version__,
        "pysam": pysam.__version__
    }
}

with open("versions.yml", "w") as f:
    yaml.dump(versions, f)
