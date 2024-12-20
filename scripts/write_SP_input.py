import os

from collections import defaultdict
from Levenshtein import distance
from itertools import groupby

import csv
import json

folder_path = './test/data/'

def read_sequences(file):
    with open(file, 'r') as file:
            fasta_dict = {}
            for key, group in groupby(file, lambda line: line.startswith('>')):
                if key:  # This is the header line
                    header = next(group).strip()[1:]  # Remove '>' and strip whitespace
                    header_id = int(header)  # Convert the ID to an integer
                else:  # This is the sequence line(s)
                    sequence = ''.join(line.strip() for line in group)  # Join all sequence lines
                    fasta_dict[header_id] = sequence
    return fasta_dict

def read_clusters(file):
    grouped_data = defaultdict(list)
    with open(file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        
        for row in reader:
            first_col, second_col = int(row[0]), int(row[1])
            grouped_data[first_col].append(second_col)

    # Convert the grouped data into a list of lists
    result = [sorted(group) for group in grouped_data.values()]
    return result

def compute_pairwise_distances(clusters, sequences):
    # Example dummy distance function (replace with real computation)
    pairwise_distances = []

    for cluster in clusters:
        cluster_distances = {}
        for i, seq1 in enumerate(cluster):
            for seq2 in cluster[i+1:]:  # Only compute when seq1 < seq2
                dist = distance(sequences[seq1], sequences[seq2])
                if seq1 not in cluster_distances:
                    cluster_distances[seq1] = {}
                cluster_distances[seq1][seq2] = dist
        pairwise_distances.append(cluster_distances)

    return pairwise_distances

def create_cluster_map(clusters):
    # Create a mapping of sequence ID -> cluster index
    cluster_map = {}
    for cluster_idx, cluster in enumerate(clusters):
        for seq_id in cluster:
            cluster_map[seq_id] = cluster_idx + 1
    return cluster_map

def save_cluster_info(out_file, clusterfile, sequencesfile):

    clusters = read_clusters(clusterfile)
    sequences = read_sequences(sequencesfile)

    data = {
        "cluster_map": create_cluster_map(clusters),
        "distances": compute_pairwise_distances(clusters, sequences)
    }
    
    with open(out_file, 'w') as f:
        json.dump(data, f, indent=4)


names = []
with open(folder_path + 'workfiles/names.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        names.append(row[0])

valid_groups = {}
with open(folder_path + 'workfiles/valid_groups.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        valid_groups[row['species']] = json.loads(row['accesions']) 

name_to_index = {name: idx for idx, name in enumerate(names)}

with open(folder_path + 'workfiles/feat_seqs.csv', 'r') as csvfile:
    reader = csv.reader(csvfile)
    feat_seqs = [row for row in reader]


os.mkdir(folder_path + 'workfiles/sp_fasta')
os.mkdir(folder_path + 'workfiles/sp_json')
for file in os.listdir(folder_path + "workfiles/groups"):
    with open(folder_path + "workfiles/groups/" + file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        groups = [row for row in reader]

    seq_dict = {}
    counter = 1
    for v in valid_groups[file[:-4]]:
        for seq in feat_seqs[name_to_index[v]]:
            if seq not in seq_dict:
                seq_dict[seq] = counter
                counter += 1

    species = file[:-4].replace(' ', '_')

    count = 0
    for group in groups:
        if len(group) > 4:

            save_cluster_info(folder_path + 'workfiles/sp_json/' + species + '_' + str(count) + '.json', folder_path + 'cluster/cluster_tsv/' + species + '.tsv', folder_path + 'cluster/cluster_fasta/' + species + '.fa')

            file = open(folder_path + "workfiles/sp_fasta/" + species + "_" + str(count) + ".fa", "w")
            for v in group:
                file.write(">" + v + "\n")
                for seq in feat_seqs[name_to_index[v]]:
                    file.write(str(seq_dict[seq]) + ", ")
                file.seek(file.tell() - 2, 0)
                file.truncate()
                file.write("\n")

            file.close()

