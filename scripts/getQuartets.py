import numpy as np 
import pandas as pd
import os
import glob
import json
from tqdm import tqdm
from pathlib import Path
from Bio import Phylo
import treeswift
from io import StringIO
import sys

def get_all_csv_paths(folder_path):
    all_paths = glob.glob(folder_path + "/*/*.csv") + glob.glob(folder_path + "/*.csv")
    return all_paths

# all_paths = get_all_csv_paths("example/simulated_data")

def get_sorted_tuple(tup):
    return str(sorted(tup))

def get_quartets(csv_path: str, 
                # use_original_weighting = True, # whether to use the "weight" field in the original CSV (false sets all weights per character to 1)
    ):
    df = pd.read_csv(csv_path)
    all_quartets = {} # string to weight
    no_quartet_count = 0 
    # for all characters
    for row in df.iterrows():
        # weight = 1.0 if not use_original_weighting else float(row[1].iloc[2]) # morphological characters weighted > lexical & phonological characters
        row = row[1][3:].to_dict()
        # row is a dictionary of language -> exhibit
        exhibits_state = {} 
        # exhibits_gene is state -> which languages have this state 
        for taxon, states in row.items():
            states = str(states).split('/')
            for c in states:
                if c not in exhibits_state:
                    exhibits_state[c] = []
                exhibits_state[c].append(taxon)
        # keys: all states with more than 2 languages exhibiting it 
        keys = [
            k 
            for (k, v) in exhibits_state.items() 
            if len(v) > 1 
        ]
        quartets = {} # quartet (canonical string) -> dict (tuple) -> count
        for i in range(len(keys)):
            for j in range(i):
                # l: the lists of the languages exhibiting states i and j
                l1 = exhibits_state[keys[i]]
                l2 = exhibits_state[keys[j]]
                l1p = [(a, b) for a in l1 for b in l1 if a < b] # a < b so that only AB gets pushed in (and not BA) 
                l2p = [(a, b) for a in l2 for b in l2 if a < b]
                all_q = [ # all possible quartets supported by this character
                    (a, b, c, d) for (a, b) in l1p for (c, d) in l2p 
                    if (a != c and a != d and b != c and b != d)
                ]
                # the conditional is so that you can't have things like AB|AC

                for q in all_q: 
                    canonical_form = get_sorted_tuple(q) # this is the canonical form of the quartet, depending only on what leaves it has
                    if canonical_form not in quartets: 
                        quartets[canonical_form] = {}
                    if q not in quartets[canonical_form]:
                        quartets[canonical_form][q] = 0
                    quartets[canonical_form][q] += 1

        # get most likely quartet for each four languages
        
        for _, four_quartets in quartets.items():
            list_form = sorted([
                (frq, q) for (q, frq) in four_quartets.items()
            ], reverse=True) # a list of (frequency of the quartet, quartet) trees sorted in decreasing order of frequency
            if len(list_form) > 1 and list_form[0][0] == list_form[1][0]:
                # throw away these four leaves if the top two counts are the same
                no_quartet_count += 1
                continue 
            returned_quartet = list_form[0][1]
            if returned_quartet not in all_quartets:
                all_quartets[returned_quartet] = 0
            all_quartets[returned_quartet] += 1
    # print(f"In {no_quartet_count} instances, two quartets were equally likely for a character & four leaves.")
    return (no_quartet_count, all_quartets)

def print_quartets(csv_path: str): 
    _, quartets = get_quartets(csv_path=csv_path)
    # print("HI!!!")
    for q, w in quartets.items():
        (a, b, c, d) = q
        print(f'(({a},{b}),({c},{d}));\n' * w, end="") 

def make_sure_path_exists(path):
    path = '/'.join(path.split('/')[:-1])
    # print("PATH = ", path)
    Path(path).mkdir(parents=True, exist_ok=True)

def get_quartets_from_folder(csvs_folder: str, save_metadata: bool, output_folder: str = "."):
    all_paths = get_all_csv_paths(csvs_folder)
    print(all_paths) 
    metadatas = []

    for path in tqdm(all_paths):
        nqc, quartets = get_quartets(path)
        input_path_ASTRAL = output_folder + "/input_ASTRAL/" + path.split(csvs_folder)[1][:-4] + ".tre"
        # input_path_wQMC = output_folder + "/input_wQMC/" + path.split(csvs_folder)[1][:-4] + ".txt"
        md = {
            "src": path,
            "input_ASTRAL": input_path_ASTRAL,
            # "input_wQMC": input_path_wQMC,
            "no_quartet_count": nqc
        }
        # print(os.dirname(output_path_ASTRAL))
        make_sure_path_exists(input_path_ASTRAL)
        # make_sure_path_exists(input_path_wQMC)
        with open(input_path_ASTRAL, "w") as f:
            for q, w in quartets.items():
                (a, b, c, d) = q
                f.write(f'(({a},{b}),({c},{d}));\n') 
        # with open(input_path_wQMC, "w") as f:
        #     for q, w in quartets.items():
        #         (a, b, c, d) = q
        #         f.write(f'{a},{b}|{c},{d}:{w}\n') 
        metadatas.append(md)
    if save_metadata:
        with open("./metadata.json", "w") as f:
            json.dump(metadatas, f)

def nexus_to_newick(trees_path): # expects a path to a .trees file under GA/ in NEXUS format. Reads all trees.
    with open(trees_path, "r") as f:
        data = Phylo.NexusIO.parse(f)
        output = StringIO()
        written_trees = Phylo.NewickIO.write(data, output, plain=True) # plain = True removes branch lengths which aren't important
        newick_trees = output.getvalue().split('\n')[:-1] # there's a newline following the last tree
        def resolve_tree(newick_tree):
            tree = treeswift.read_tree(newick_tree, schema='newick')
            #print(newick_tree)
            tree.resolve_polytomies()
            return tree.newick()
        trees_resolved = list(map(
        resolve_tree 
        , newick_trees))
        output.close()
        for t in trees_resolved:
            print(t)
        print(f'{len(trees_resolved)} Trees Written with arbitrary polytomy splitting', file=sys.stderr)