import argparse

from scripts.functions import *


def main():
    parser = argparse.ArgumentParser(description="Run the genome processing pipeline.")
    parser.add_argument(
        "input_dir",
        help="Path to work directory (same path as used for fetch_data.py). (Positional argument)"
    )
    args = parser.parse_args()

    # Paths from arguments
    input_dir = args.input_dir
    workfolder = os.path.join(input_dir, 'workfolder')
    os.mkdir(workfolder)

    # extract trees
    extract_trees(workfolder)

    # fix trees
    fix_trees(workfolder)


if __name__ == "__main__":
    main()
