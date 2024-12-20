import argparse
import subprocess

from scripts.functions import *


def main():
    parser = argparse.ArgumentParser(description="Run the genome processing pipeline.")
    parser.add_argument(
        "input_dir",
        help="Path to directory containing the input directory with the .txt files containing the RefSeq accessions. (Positional argument)"
    )
    parser.add_argument('-t', '--trees', action='store_true',
                        help='Download gbk files for panX')
    args = parser.parse_args()

    # Paths from arguments
    input_dir = args.input_dir
    trees = args.trees
    workfolder = os.path.join(input_dir, 'workfolder')
    os.mkdir(workfolder)

    # download genome fasta files
    try:
        subprocess.run(
            ["bash", "./scripts/download_ncbi_fasta.sh", "-i", input_dir + '/input', "-o", workfolder],
            check=True,
        )
        print("Bash script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing script: {e}")

    # run IntegronFinder
    try:
        subprocess.run(
            ["bash", "./scripts/run_IF.sh", "-p", workfolder],
            check=True,
        )
        print("Bash script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing script: {e}")

    # extract IF output
    extract_IF_output(workfolder)

    # extract integron data and write csv
    feat_seqs, names, valid_groups = read_IF_write_csv(workfolder)

    # run cluster
    try:
        subprocess.run(
            ["bash", "./scripts/run_cluster.sh", "-d", os.path.join(workfolder, 'cluster', 'cluster_fasta')],
            check=True,
        )
        print("Bash script executed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing script: {e}")

    # overlap groups
    overlap_groups(workfolder, feat_seqs, names, valid_groups)

    # download gbk files for panX
    if trees:
        try:
            subprocess.run(
                ["bash", "./scripts/download_ncbi_to_panX.sh", "-i", os.path.join(workfolder, 'panX_txt'), "-o",
                 os.path.join(workfolder, 'panX')],
                check=True,
            )
            print("Bash script executed successfully.")
        except subprocess.CalledProcessError as e:
            print(f"Error executing script: {e}")

    # write SpacerPlacer input files
    write_sp_files(workfolder, feat_seqs, names, valid_groups)


if __name__ == "__main__":
    main()
