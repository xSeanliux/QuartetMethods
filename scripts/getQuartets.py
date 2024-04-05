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
from functools import reduce
from .loss_quartets import *
from itertools import combinations

def get_all_csv_paths(folder_path):
    all_paths = glob.glob(folder_path + "/*/*.csv") + glob.glob(folder_path + "/*.csv")
    return all_paths

# all_paths = get_all_csv_paths("example/simulated_data")

def get_sorted_tuple(tup):
    return str(sorted(tup))

class Frequenter():
    # supports initialisation with a list of items
    # keeps a frequency table 
    # supports querying of a list of items, returning the most frquent out of them. 
    # ties are broken arbitrarily.
    def __init__(self, arr):
        self.frqs = {}
        for x in arr:
            if x not in self.frqs:
                self.frqs[x] = 0
            self.frqs[x] += 1
    
    def query(self, k):
        return max(k, key=lambda x : self.frqs[x])

def process_row(row):
    data = row[3:]
    data = [x.split('/') for x in data]
    frq = Frequenter(reduce(
        lambda acc, x : acc + x,
        data,
        [] 
    ))
    data = [ frq.query(k) for k in data ]
    row[3:] = data
    return row

def resolve_polymorphism_using_mp4(df):
    for i, row in enumerate(df.values):
        df.iloc[i, :] = process_row(row)
    return df

def get_loss_based_quartets(df, mode):
    assert(4 <= mode <= 5)
    quartets = {}
    names = df.columns.values[3:]
    values = df.values[:, 3:]
    for row in values: 
        states_dict = {
            n: s.split('/') for (n, s) in zip(names, row)
        }
        #print(states_dict)
        for nametup in combinations(names, 4):
            best_trees = get_best_quartets_loss({
                x : states_dict[x] for x in nametup
            })

            if(len(best_trees) >= 15): # when all trees have the same score, no information is supplied.
                continue
            
            # list of best trees -> get quartets from them -> turn back into a list
            if(mode == 4):
                best_trees = list(map(lambda x: x.rooted_quartet(), best_trees))
            elif(mode == 5):
                best_trees = list(map(lambda x: x.get_unrooted_tuple(), best_trees))

            for tup in best_trees:
                if tup not in quartets:
                    quartets[tup] = 0
                quartets[tup] += 1
    return (None, quartets) 
    # return (None, quartets)
    # where quartets is a dict (a,b,c,d) -> w

def get_quartets(csv_path   : str, 
                mode        : int=1
    ):
    '''
    Mode 
    1 Only most probable
    2 All most probable if they're not equally probable (don't do this one)
    3 resolve polymorphism using the MP4 method (most frequent appearing exhibited state)
    4 Loss-based, rooted (with extra root leaf)
    5 Loss-based, unrooted (just give quartets)
    '''
    df = pd.read_csv(csv_path)
    if mode == 3:
        df = resolve_polymorphism_using_mp4(df)
    if 4 <= mode <= 5:
        return get_loss_based_quartets(df=df, mode=mode)

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
                l1p = list(combinations(l1, 2))
                l2p = list(combinations(l2, 2)) 
                all_q = [ # all possible quartets supported by this character
                    (a, b, c, d) for (a, b) in l1p for (c, d) in l2p 
                    if (a != c and a != d and b != c and b != d)
                ]
                # the conditional is so that you can't have things like AB|AC
                # kill `$(ps aux | grep 'zxliu2/.' | awk '{print $2}')`
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
            most_common_cnt = list_form[0][0]
            most_common_quartets = list(filter(lambda x: x[0] == most_common_cnt, list_form)) # get just the most common quartets
            if(mode == 1 or mode == 3):
                if len(most_common_quartets) > 1:
                    continue
                returned_quartet = most_common_quartets[0][1]
                if returned_quartet not in all_quartets:
                    all_quartets[returned_quartet] = 0
                all_quartets[returned_quartet] += 1
            elif(mode == 2): # deprecated but this was implemented already so it stays
                if len(most_common_quartets) >= 3:
                    continue
                for _, q in most_common_quartets:
                    returned_quartet = q
                    if returned_quartet not in all_quartets:
                        all_quartets[returned_quartet] = 0
                    all_quartets[returned_quartet] += 1

                
    # print(f"In {no_quartet_count} instances, two quartets were equally likely for a character & four leaves.")
    return (no_quartet_count, all_quartets)

def print_quartets(csv_path: str, mode: int): 
    _, quartets = get_quartets(csv_path=csv_path, mode=mode)
    # print("HI!!!")
    for q, w in quartets.items():
        if(type(q) == tuple):
            (a, b, c, d) = q
            print(f'(({a},{b}),({c},{d}));\n' * w, end="") 
        elif(type(q) == str):
            if(q[-1] != '\n'):
                q = q + '\n'
            print(q * w, end="")

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
