#!/bin/bash

usage() {
    echo "Usage: $0 -d DIRECTORY"
    exit 1
}

while getopts "d:" opt; do
    case "$opt" in
        d) INPUT_DIR=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" ]]; then
    usage
fi

if [[ ! -d "$INPUT_DIR" ]]; then
    echo "Error: Directory $INPUT_DIR does not exist."
    exit 1
fi

for file in "$INPUT_DIR"/*.fa; do
    base=$(basename "$file" .fa)

    echo "Processing $base"

    base_dir="$(dirname "$INPUT_DIR")/$base"
    mkdir -p "$base_dir"

    mmseqs createdb "$file" "$base_dir/$base.db"

    mmseqs cluster "$base_dir/$base.db" "$base_dir/$base.clu" "$base_dir/tmp" --min-seq-id 0.8 -c 0.8

    mmseqs createtsv "$base_dir/$base.db" "$base_dir/$base.db" "$base_dir/$base.clu" "$base_dir/$base.tsv"

    rm -rf "$base_dir/tmp"
done
