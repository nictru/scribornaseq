#!/usr/bin/env python

# Set numba chache dir to current working directory (which is a writable mount also in containers)
import os

os.environ["NUMBA_CACHE_DIR"] = "."

import scanpy as sc
import pandas as pd
import platform
from scipy.sparse import csr_matrix

def format_yaml_like(data: dict, indent: int = 0) -> str:
    """Formats a dictionary to a YAML-like string.
    Args:
        data (dict): The dictionary to format.
        indent (int): The current indentation level.
    Returns:
        str: A string formatted as YAML.
    """
    yaml_str = ""
    for key, value in data.items():
        spaces = "  " * indent
        if isinstance(value, dict):
            yaml_str += f"{spaces}{key}:\\n{format_yaml_like(value, indent + 1)}"
        else:
            yaml_str += f"{spaces}{key}: {value}\\n"
    return yaml_str

def dump_versions():
    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "scanpy": sc.__version__,
            "pandas": pd.__version__
        }
    }

    with open("versions.yml", "w") as f:
        f.write(format_yaml_like(versions))

#
# Run main script
#

adata = sc.read_mtx("${input}/matrix.mtx.gz").T
adata.obs = pd.read_csv("${input}/barcodes.tsv.gz", header=None, sep="\\t")
adata.var = pd.read_csv("${input}/features.tsv.gz", header=None, sep="\\t")
adata.var.columns = ["chr", "start", "end", "strand", "motif", "annotated", "unique_reads", "multi_reads", "max_overhang"]
adata.var["strand"] = adata.var["strand"].map({0: "undefined", 1: "+", 2: "-"})
adata.var["motif"] = adata.var["motif"].map({0: "non-canonical", 1: "GT/AG", 2: "CT/AC", 3: "GC/AG", 4: "GT/GC", 5: "AT/AC", 6: "GT/AT"})
adata.var["annotated"] = adata.var["annotated"].map({0: "No", 1: "Yes"})
adata.var.index = adata.var["chr"] + ":" + adata.var["start"].astype(str) + "-" + adata.var["end"].astype(str) + ":" + adata.var["strand"]
adata.obs.columns = ["barcode"]
adata.X = csr_matrix(adata.X)
adata.obs["sample"] = "${meta.id}"

adata.write_h5ad("${meta.id}_matrix.h5ad")

# dump versions
dump_versions()
