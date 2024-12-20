#!/bin/bash

usage() {
    echo "Usage: $0 -d DATA_DIR -p PANX_PATH"
    exit 1
}

while getopts "d:p:e:" opt; do
    case "$opt" in
        d) DATA_DIR=$OPTARG ;;
        p) PANX_PATH=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ -z "$DATA_DIR" || -z "$PANX_PATH" ]]; then
    usage
fi

if [[ ! -d "$DATA_DIR" ]]; then
    echo "Error: Directory $DATA_DIR does not exist."
    exit 1
fi

if [[ ! -f "$PANX_PATH" ]]; then
    echo "Error: panX.py script not found at $PANX_PATH."
    exit 1
fi

for subfolder in "$DATA_DIR"/*/; do
    subfolder_name=$(basename "$subfolder")

    echo "Activating Conda environment $ENV_NAME..."
    (
        echo "Running $PANX_PATH for subfolder $subfolder_name..." && \
        python "$PANX_PATH" -fn "$subfolder" -sl "$subfolder_name" -t 4
    )

    if [[ $? -eq 0 ]]; then
        echo "Successfully processed $subfolder_name"
    else
        echo "Error processing $subfolder_name"
    fi
done

rm -rf tmp

echo "All subfolders processed."
