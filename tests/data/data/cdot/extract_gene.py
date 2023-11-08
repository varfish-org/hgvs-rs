#!/usr/bin/env python
"""Extract transcripts from cdot JSON file."""

import gzip
import json
import os
import sys
import typing

#: HGNC identifier ofs gene to extract, from env or BRCA1.
HGNC_IDS = os.environ.get("HGNC_IDS", "HGNC:1100").split(",")


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
    print("gettings keys...", file=sys.stderr)
    keys = []
    for key, gene in full_data["genes"].items():
        if "hgnc" in gene and gene["hgnc"] in hgncs_stripped and key not in keys:
            keys.append(key)
    print(f"- keys = {keys}", file=sys.stderr)
    print("extracting...", file=sys.stderr)
    full_data["genes"] = {key: full_data["genes"][key] for key in keys}
    full_data["transcripts"] = {
        key: value
        for key, value in full_data["transcripts"].items()
        if "hgnc" in value and value["hgnc"] in hgncs_stripped
    }
    return full_data


def main(json_paths: typing.List[str]):
    """Extract transcripts from cdot JSON file."""
    hgncs_stripped = [hgnc_id.replace("HGNC:", "") for hgnc_id in HGNC_IDS]
    for json_path in json_paths:
        gene_data = extract_gene(json_path, hgncs_stripped)
        symbols = "_".join(
            [gene["gene_symbol"].lower() for gene in gene_data["genes"].values()]
        )
        out_path = json_path.replace(".gz", "").replace(".json", f".{symbols}.json")
        print(f"writing to {out_path}...", file=sys.stderr)
        with open(out_path, "wt") as outputf:
            json.dump(gene_data, outputf, indent=2)


if __name__ == "__main__":
    main(sys.argv[1:])
