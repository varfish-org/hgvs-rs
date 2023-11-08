#!/usr/bin/env python
"""Extract transcripts from cdot JSON file."""

import gzip
import json
import os
import sys
import typing

#: HGNC identifier ofs gene to extract, from env or BRCA1.
HGNC_IDS = os.environ.get("HGNC_ID", "HGNC:672").split(",")


def load_json(json_path: str) -> typing.Dict[str, typing.Any]:
    """Load JSON file."""
    print(f"Loading {json_path}", file=sys.stderr)
    if json_path.endswith(".gz"):
        with gzip.open(json_path, "rt") as json_file:
            return json.load(json_file)
    else:
        with open(json_path, "rt") as json_file:
            return json.load(json_file)


def extract_gene(
    json_path: str, hgncs_stripped: typing.List[str]
) -> typing.Dict[str, typing.Any]:
    """Extract data of one gene, specified by HGNC identifier without prefix."""
    full_data = load_json(json_path)
    print("extracting...", file=sys.stderr)
    full_data["genes"] = {key: full_data["genes"][key] for key in hgncs_stripped}
    gene_symbols = [
        full_data["genes"][hgnc_id]["gene_symbol"] for hgnc_id in hgncs_stripped
    ]
    full_data["transcripts"] = {
        key: value
        for key, value in full_data["transcripts"].items()
        if value["gene_name"] in gene_symbols
    }
    return full_data


def main(json_paths: typing.List[str]):
    """Extract transcripts from cdot JSON file."""
    hgncs_stripped = [hgnc_id.replace("HGNC:", "") for hgnc_id in HGNC_IDS]
    for json_path in json_paths:
        gene_data = extract_gene(json_path, hgncs_stripped)
        symbols = [
            gene_data["genes"][hgnc_id]["gene_symbol"].lower()
            for hgnc_id in hgncs_stripped
        ].join("_")
        out_path = json_path.replace(".gz", "").replace(".json", f".{symbols}.json")
        print(f"writing to {out_path}...", file=sys.stderr)
        with open(out_path, "wt") as outputf:
            json.dump(gene_data, outputf, indent=2)


if __name__ == "__main__":
    main(sys.argv[1:])
