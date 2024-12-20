#!/bin/bash

INPUT_DIR=""
OUTPUT_DIR=""

usage() {
    echo "Usage: $0 -i INPUT_DIR -o OUTPUT_DIR"
    exit 1
}

while getopts "i:o:" opt; do
    case "$opt" in
        i) INPUT_DIR=$OPTARG ;;
        o) OUTPUT_DIR=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    usage
fi

for file in "$INPUT_DIR"/*.txt; do
    base_name=$(basename "$file" .txt)
    base_folder="$OUTPUT_DIR/$base_name"

    datasets download genome accession --inputfile "$file" --include genome --assembly-level complete

    if [ -f ncbi_dataset.zip ]; then
        mkdir -p "$base_folder"
        unzip -q ncbi_dataset.zip -d "$base_folder"
        rm ncbi_dataset.zip
    else
        echo "Error: ncbi_dataset.zip not found for $file"
        continue
    fi

    for folder in "$base_folder"/ncbi_dataset/data/*/; do
        if [ -d "$folder" ]; then
            folder_name=$(basename "$folder")
            file=$(find "$folder" -maxdepth 1 -type f | head -n 1)
            if [ -n "$file" ]; then
                extension="${file##*.}"
                mv "$file" "$folder/${folder_name}.${extension}"
            fi
        fi
    done

    mkdir -p "$base_folder/genomes"
    find "$base_folder"/ncbi_dataset/data/ -name "*.fna" -exec cp {} "$base_folder/genomes/" \;

    find "$base_folder" -mindepth 1 -maxdepth 1 ! -name "genomes" -exec rm -rf {} +
done
