import csv
import json
import os
import shutil
from collections import defaultdict
from collections import Counter
from Levenshtein import distance
from itertools import groupby
from Bio import Phylo
from Bio import SeqIO


def extract_IF_output(folder_path):
    strain_counts = defaultdict(int)
    strain_appeared = defaultdict(bool)
    for species in os.listdir(folder_path):
        os.mkdir(os.path.join(folder_path, species, 'Integrons'))
        for strain in os.listdir(os.path.join(folder_path, species, "IF_out")):
            subdir_path = os.path.join(folder_path, species, "IF_out", strain)

            for filename in os.listdir(subdir_path):
                file_path = os.path.join(subdir_path, filename)
                if not strain_appeared[strain]:
                    new_filename = f"{strain}{os.path.splitext(filename)[1]}"
                    strain_appeared[strain] = True
                else:
                    strain_counts[strain] += 1

                    new_filename = f"{strain}_{strain_counts[strain]}{os.path.splitext(filename)[1]}"

                new_file_path = os.path.join(folder_path, species, 'Integrons', new_filename)
                shutil.copy(file_path, new_file_path)


def read_IF_write_csv(folder_path):
    features = []
    names = []
    DNA_seqs = []
    names_count = {}
    reverse = []
    id_to_name = {}

    for species in os.listdir(folder_path):
        integrons_path = os.path.join(folder_path, species, 'Integrons')
        file_list = sorted(os.listdir(integrons_path))
        for filename in file_list:
            if filename.endswith(".gbk"):

                x = [rec for rec in SeqIO.parse(os.path.join(integrons_path, filename), "genbank")]

                names_count[filename[:-4]] = 0
                id_to_name[filename[:-4]] = species

                complete = False
                feats = []

                temp_feats = []

                must_reverse = False

                for rec in x:
                    for feat in rec.features:
                        if feat.type == "integron":
                            if len(temp_feats) > 0:
                                feats.append(temp_feats)
                                temp_feats = []

                            if len(feats) > 0:
                                features.append(feats)
                                DNA_seqs.append(rec.seq)
                                names.append(filename[:-4] + "_" + str(names_count[filename[:-4]]))
                                names_count[filename[:-4]] += 1
                                feats = []
                                reverse.append(must_reverse)

                            if "complete" in feat.qualifiers["integron_type"]:
                                complete = True
                            else:
                                complete = False

                        if complete and feat.type == "integrase":
                            if len(temp_feats) > 0 or len(feats) > 0:
                                must_reverse = True
                            else:
                                must_reverse = False

                        if complete and feat.type == "CDS":
                            temp_feats.append(feat)

                        if complete and feat.type == "attC" and len(temp_feats) > 0:
                            feats.append(temp_feats)
                            temp_feats = []

                if len(temp_feats) > 0:
                    feats.append(temp_feats)

                if len(feats) > 0:
                    features.append(feats)
                    DNA_seqs.append(rec.seq)
                    names.append(filename[:-4] + "_" + str(names_count[filename[:-4]]))
                    reverse.append(must_reverse)

    for i in range(len(features)):
        if reverse[i]:
            features[i] = features[i][::-1]

    feat_seqs = []

    block_lengths = {}

    for i in range(len(DNA_seqs)):
        temp = []

        for feat in features[i]:

            if len(feat) not in block_lengths:
                block_lengths[len(feat)] = 1
            else:
                block_lengths[len(feat)] += 1

            block = ""
            for f in feat:
                loc = f.location
                start, end = loc.start, loc.end
                block = block + (f.qualifiers["translation"][0])[:-1]
            temp.append(block)
        feat_seqs.append(temp)

    groups = {}

    for name in names:
        if name[:15] not in id_to_name:
            print(name[:15])
            continue
        group = id_to_name[name[:15]]
        if group not in groups:
            groups[group] = [name]
        else:
            groups[group].append(name)

    valid_groups = {}

    for key, val in groups.items():
        if len(val) > 5:
            valid_groups[key] = val

    list_to_names = {}
    for lst, name in zip(feat_seqs, names):
        lst_tuple = tuple(lst)

        if lst_tuple in list_to_names:
            list_to_names[lst_tuple].append(name)
        else:
            list_to_names[lst_tuple] = [name]

    unique_arrs = {}
    for k, v in list_to_names.items():
        unique_arrs[v[0]] = k

    new_valids = {}
    for k, v in valid_groups.items():
        new_valids[k] = []
        for i in v:
            if i in unique_arrs:
                new_valids[k].append(i)

    valid_groups = {}

    for key, val in new_valids.items():
        if len(val) > 5:
            valid_groups[key] = val

    with open(os.path.join(folder_path, 'feat_seqs.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(feat_seqs)

    with open(os.path.join(folder_path, 'names.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        for string in names:
            writer.writerow([string])

    with open(os.path.join(folder_path, 'valid_groups.csv'), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['species', 'accesions'])
        for key, value in valid_groups.items():
            writer.writerow([key, json.dumps(value)])

    os.makedirs(os.path.join(folder_path, 'cluster', 'cluster_fasta'))
    name_to_index = {name: idx for idx, name in enumerate(names)}

    for key, val in valid_groups.items():

        key = key.replace(" ", "_")

        seq_dict = {}

        counter = 1

        for v in val:
            for seq in feat_seqs[name_to_index[v]]:
                if seq not in seq_dict:
                    seq_dict[seq] = counter
                    counter += 1

        file = open(os.path.join(folder_path, 'cluster', 'cluster_fasta', key + ".fa"), "w")

        counter = 0

        for key1, val1 in seq_dict.items():
            file.write(">" + str(val1) + "\n")
            file.write(key1)
            file.write("\n")

        file.close()

    return feat_seqs, names, valid_groups


def overlap_groups(folder_path, feat_seqs, names, valid_groups):
    os.mkdir(os.path.join(folder_path, 'cluster', 'cluster_tsv'))
    for dirpath, dirnames, filenames in os.walk(os.path.join(folder_path, 'cluster')):
        for filename in filenames:
            if filename.endswith(".tsv"):
                file_path = os.path.join(dirpath, filename)
                new_file_path = os.path.join(folder_path, 'cluster', 'cluster_tsv', filename)
                if not file_path == new_file_path:
                    shutil.copy(file_path, new_file_path)

    name_to_index = {name: idx for idx, name in enumerate(names)}

    os.mkdir(os.path.join(folder_path, 'groups'))
    os.mkdir(os.path.join(folder_path, 'panX_txt'))
    for file in os.listdir(os.path.join(folder_path, 'cluster', 'cluster_tsv')):
        grouped_data = defaultdict(list)

        with open(os.path.join(folder_path, 'cluster', 'cluster_tsv', file), 'r') as tsvfile:
            reader = csv.reader(tsvfile, delimiter='\t')

            for row in reader:
                first_col, second_col = int(row[0]), int(row[1])
                grouped_data[first_col].append(second_col)

        result = [group for group in grouped_data.values()]

        cluster_dict = {}
        for cluster in result:
            for seq in cluster:
                cluster_dict[seq] = result.index(cluster) + 1

        species = file[:-4]

        if species not in valid_groups:
            continue

        vals = valid_groups[species]

        seq_dict = {}
        counter = 1
        for v in vals:
            for seq in feat_seqs[name_to_index[v]]:
                if seq not in seq_dict:
                    seq_dict[seq] = cluster_dict[counter]
                    counter += 1

        numbers = []
        for v in vals:
            temp = []
            for seq in feat_seqs[name_to_index[v]]:
                temp.append(seq_dict[seq])
            numbers.append(set(temp))

        clusters = [[val] for val in vals]

        change = True

        while change:
            change = False
            for i in range(len(clusters)):
                for j in range(i + 1, len(clusters)):
                    if numbers[i].intersection(numbers[j]):
                        clusters[i].extend(clusters[j])
                        clusters.pop(j)
                        numbers[i] = numbers[i].union(numbers[j])
                        numbers.pop(j)
                        change = True
                        break
                    if change:
                        break
                if change:
                    break

        clusters = [cluster for cluster in clusters if len(cluster) > 5]

        with open(os.path.join(folder_path, 'groups', species + '.csv'), 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerows(clusters)

        count = 0
        for cluster in clusters:
            group = [acc[:15] for acc in cluster]
            group = list(set(group))
            with open(os.path.join(folder_path, 'panX_txt', species + '_' + str(count) + '.txt'), 'a') as panX_file:
                for item in group:
                    panX_file.write("%s\n" % item)
            count += 1


def extract_trees(folder_path):
    os.mkdir(os.path.join(folder_path, 'trees'))

    for folder in os.listdir(os.path.join(folder_path, 'panX')):
        if folder.endswith(".zip") or folder.endswith(".sh") or folder == "tmp":
            continue
        if os.path.exists(os.path.join(folder_path, 'panX', folder, 'vis', 'geneCluster', "strain_tree.nwk")):
            file = os.path.join(folder_path, 'panX', folder, 'vis', 'geneCluster', "strain_tree.nwk")
        else:
            file = os.path.join(folder_path, 'panX', folder, 'geneCluster', "strain_tree.nwk")

        shutil.copy(file, os.path.join(folder_path, 'trees', folder + ".nwk"))


def duplicate_in_binary_chain(parent, node, count):
    og_length = node.branch_length
    dist = og_length / (count + 1)
    parent.clades.remove(node)

    last_node = node.__class__(name=node.name, branch_length=dist)
    for i in range(count):
        new_internal_clade = node.__class__(branch_length=dist)

        duplicate_node = node.__class__(name=f"{node.name}dup{i + 1}", branch_length=dist * (i + 1))

        new_internal_clade.clades.append(last_node)
        new_internal_clade.clades.append(duplicate_node)

        last_node = new_internal_clade

    parent.clades.append(last_node)


def fix_trees(folder_path):
    for file in os.listdir(folder_path + 'groups'):
        with open(folder_path + 'groups/' + file, 'r') as csvfile:
            reader = csv.reader(csvfile)
            groups = [row for row in reader]

            groups = [group for group in groups if len(group) > 4]

            counter = 0

            for group in groups:
                group_mod = [acc[:15] for acc in group]

                tree_path = folder_path + 'trees/' + file[:-4].replace(" ", "_") + "_" + str(counter) + ".nwk"
                tree = Phylo.read(tree_path, "newick")

                if len(group_mod) != len(set(group_mod)):
                    counts = Counter(group_mod)
                    duplication_dict = {key: value for key, value in counts.items() if value > 1}

                    for clade in tree.find_clades(order="level"):
                        if clade.name in duplication_dict:
                            for parent in tree.get_nonterminals():
                                if clade in parent.clades:
                                    duplicate_in_binary_chain(parent, clade, duplication_dict[clade.name] - 1)
                                    break

                for acc in group:
                    for clade in tree.find_clades():
                        if clade.name:
                            if clade.name.count("_") <= 1:
                                if clade.name[:15] in acc:
                                    clade.name = acc
                                    break

                with open(tree_path, 'w') as output_tree_file:
                    Phylo.write(tree, output_tree_file, "newick")

                counter += 1


def read_sequences(file):
    with open(file, 'r') as file:
        fasta_dict = {}
        for key, group in groupby(file, lambda line: line.startswith('>')):
            if key:
                header = next(group).strip()[1:]
                header_id = int(header)
            else:
                sequence = ''.join(line.strip() for line in group)
                fasta_dict[header_id] = sequence
    return fasta_dict


def read_clusters(file):
    grouped_data = defaultdict(list)
    with open(file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')

        for row in reader:
            first_col, second_col = int(row[0]), int(row[1])
            grouped_data[first_col].append(second_col)

    result = [sorted(group) for group in grouped_data.values()]
    return result


def compute_pairwise_distances(clusters, sequences):
    pairwise_distances = []

    for cluster in clusters:
        cluster_distances = {}
        for i, seq1 in enumerate(cluster):
            for seq2 in cluster[i + 1:]:
                dist = distance(sequences[seq1], sequences[seq2])
                if seq1 not in cluster_distances:
                    cluster_distances[seq1] = {}
                cluster_distances[seq1][seq2] = dist
        pairwise_distances.append(cluster_distances)

    return pairwise_distances


def create_cluster_map(clusters):
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


def write_sp_files(folder_path, feat_seqs, names, valid_groups):
    name_to_index = {name: idx for idx, name in enumerate(names)}
    os.mkdir(os.path.join(folder_path, 'sp_fasta'))
    os.mkdir(os.path.join(folder_path, 'sp_json'))
    for file in os.listdir(folder_path + "/groups"):
        with open(folder_path + "/groups/" + file, 'r') as csvfile:
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

                save_cluster_info(folder_path + '/sp_json/' + species + '_' + str(count) + '.json',
                                  folder_path + '/cluster/cluster_tsv/' + species + '.tsv',
                                  folder_path + '/cluster/cluster_fasta/' + species + '.fa')

                file = open(folder_path + "/sp_fasta/" + species + "_" + str(count) + ".fa", "w")
                for v in group:
                    file.write(">" + v + "\n")
                    for seq in feat_seqs[name_to_index[v]]:
                        file.write(str(seq_dict[seq]) + ", ")
                    file.seek(file.tell() - 2, 0)
                    file.truncate()
                    file.write("\n")

                file.close()
