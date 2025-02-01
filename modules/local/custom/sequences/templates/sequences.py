#!/usr/bin/env python3

import Bio
from Bio import SeqIO
import yaml
import platform

with SeqIO.parse("${fasta}", "fasta") as records, \
    open("${prefix}.fasta", "w") as output, \
    open("${prefix}.gtf", "w") as gtf:
    for record in records:
        gene_type = record.description.split(" ")[1]

        padding_length = int(1e3)

        start = 1 + padding_length
        end = start + len(record.seq) - 1
        chr = f"cusom_{record.id}"
        identifier = record.id

        for feature_type in ["gene", "transcript", "exon"]:
            prefix = f"{chr}\\tCUSTOM\\t{feature_type}\\t{start}\\t{end}\\t.\\t+\\t.\\t"
            if feature_type == "gene":
                gtf.write(prefix + f'gene_id "{identifier}"; gene_type "{gene_type}"; gene_name "{identifier}";\\n')
            else:
                gtf.write(prefix + f'gene_id "{identifier}"; transcript_id "{identifier}.1"; gene_type "{gene_type}"; gene_name "{identifier}";\\n')

        record.id = chr
        record.description = ""
        padding = "N" * padding_length
        record.seq = padding + record.seq + padding

        SeqIO.write(record, output, "fasta")

versions = {
    "versions": {
        "biopython": Bio.__version__,
        "python": platform.python_version(),
    }
}

with open("versions.yml", "w") as output:
    yaml.dump(versions, output)
