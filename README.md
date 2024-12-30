# IntegronEvo

## Installation and Usage
IntegronEvo is supported on Linux, there is no guarantee it will work on other operating systems.  
To install and run it follow these steps:
1. Clone the GitHub repository:

    ```bash
    git clone https://github.com/Fabian1567/IntegronEvo.git
    ```

2. Navigate to the cloned repository and install requirements in a new conda environment:

    ```bash
    cd IntegronEvo
    conda env create -n <myenv> -f environment/environment.yml
    ```
    With -n <myenv> you can name the environment. You can choose any name you like or leave it out, as environment.yml provides a default name (integronevo_env).

3. Activate the environment:

   ```bash
   conda activate <myenv>
   ```

4. Create necessary SpacerPlacer input data by running:
    ```bash
   python fetch_data.py <input_path> [--trees]
   ```
   Information on how to structure input files can be found below.  
    Only use the `--trees` flag if you intend to create trees using panX instead of relying on the trees created by SpacerPlacer.

5. (Optional) If you want to create trees using panX install it as      explained [here](https://github.com/neherlab/pan-genome-analysis).  
Activate the panX conda environment and run panX on all subfolders in `input_path/workfolder/panX` or alternatively run:
    ```bash
   ./scripts/run_panX.sh -d input_path/workfolder/panX -p <path to panX.py>
   ```
   The trees created by panX require postprocessing, therefore switch back to the `<myenv>` environment and run:
   ```bash
   python fix_trees.py <input_path>
   ```

6. Install SpacerPlacer as explained [here](https://github.com/Fabian1567/SpacerPlacer_IntegronEvo) activate the SpacerPlacer conda environment and run:
    ```bash
   python <path to spacerplacer.py> input_path/workfolder/sp_fasta/ <output_path> --cluster_json input_path/workfolder/sp_json/ [--tree_path input_path/workfolder/trees/]
   ```
   Only use the `--tree_path` flag if you performed step 5!


## Input
Input has to be structured as follows:
```bash
.
└── data
    └── input
        ├── group_1.txt
        └── group_2.txt

```
Where data is the folder referred to as <input_path> above and each `.txt` file contains a list of RefSeq accesions like:
```
GCF_013487985.1
GCF_016598815.1
GCF_014158455.1
GCF_016861425.1
GCF_007179295.1
GCF_003925855.2
```
Make sure the filenames of the `.txt` files contain no special characters and spaces are replaced by underscores.
