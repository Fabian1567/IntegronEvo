#!/bin/bash

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

mkdir -p "$OUTPUT_DIR"

for file in "$INPUT_DIR"/*.txt; do
    base_name=$(basename "$file" .txt)

    datasets download genome accession --inputfile "$file" --include gbff --assembly-level complete

    if [ -f ncbi_dataset.zip ]; then
        echo "Processing $file..."

        base_output_dir="$OUTPUT_DIR/$base_name"
        mkdir -p "$base_output_dir"

        unzip ncbi_dataset.zip -d "$base_output_dir"

        rm ncbi_dataset.zip

        DATA_DIR="$base_output_dir/ncbi_dataset/data"

        if [[ -d "$DATA_DIR" ]]; then
            echo "Organizing data for $base_name..."

            INPUT_GENBANK="$base_output_dir/input_GenBank"
            mkdir -p "$INPUT_GENBANK"

            for subfolder in "$DATA_DIR"/*/; do
                folder_name=$(basename "$subfolder")

                gbff_file=$(find "$subfolder" -name "*.gbff" -type f)

                if [[ -f "$gbff_file" ]]; then
                    new_filename="$folder_name.gbk"
                    cp "$gbff_file" "$INPUT_GENBANK/$new_filename"
                else
                    echo "No .gbff file found in $subfolder"
                fi
            done

            find "$base_output_dir" -mindepth 1 -maxdepth 1 ! -name "input_GenBank" -exec rm -rf {} +
        else
            echo "No data directory found in $base_output_dir"
        fi
    else
        echo "Error: ncbi_dataset.zip not found for $file"
    fi
done

echo "All tasks completed. Each dataset now contains its own input_GenBank folder."
