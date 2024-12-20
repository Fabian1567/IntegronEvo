#!/bin/bash

PARENT_DIR=""

usage() {
    echo "Usage: $0 -p PARENT_DIR"
    exit 1
}

while getopts "p:" opt; do
    case "$opt" in
        p) PARENT_DIR=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$PARENT_DIR" ]]; then
    usage
fi

for folder in "$PARENT_DIR"/*/; do
    genomes_folder="${folder}genomes"

    if [ -d "$genomes_folder" ]; then
        if_out_folder="${folder}IF_out"
        mkdir -p "$if_out_folder"

        for file in "$genomes_folder"/*; do
            abs_file_path=$(realpath "$file")

            (
                cd "$if_out_folder" || exit
                integron_finder --local-max --func-annot --gbk --cpu 4 "$abs_file_path"
            )
        done

        (
            cd "$if_out_folder" || exit

            rename 's/Results_Integron_Finder_//' Results_Integron_Finder_*
            find . -mindepth 1 -type f ! -name "*.gbk" -delete
            find . -type d -empty -delete

            find . -type f -exec sh -c '
                for file do
                    count=$(grep -o -i "complete" "$file" | wc -l)
                    if [ "$count" -lt 2 ]; then
                        rm "$file"
                    fi
                done
            ' sh {} +
        )
    else
        echo "No genomes folder found in $folder"
    fi
done
