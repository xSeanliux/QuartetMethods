import numpy as np  

class RootedQuartet:
    '''
    This class is used to store a particular rooted binary tree on four leaves and to calculate the Dollo loss on it.
    It is more wordy but should be more efficient than actually calling Dollo, since it is specifically tailored to handle quartets.
    '''
    def __init__(self, type, labels, pair):
        # pair: a pair of two names
        '''
        Initializes the RootedQuartet.
        @params
        type: BAL and LOP (balanced, lopsided)
            BAL: ((a, b), (c, d)), described by an UNORDERED pair (a, b). 
            LOP: (a, (b, (c, d))), described by an ORDERED pair (a, b).
        labels: list of four names of the leaves
        pair: a pair of names which must be in labels, described above
        
        leaves are stored numerically as integers [0-3], self.labels is a mapping from integer to string. 
        For example, calling
        RootedQuartet("BAL", ["a", "b", "c", "d"], ("a", "b")) -> ((a, b), (c, d))
        RootedQuartet("LOP", ["a", "b", "c", "d"], ("c", "b")) -> (c, (b, (a, d)))
        '''
        assert(type in ['BAL', 'LOP'])
        assert(pair[0] in labels and pair[1] in labels)
        self.labels = labels 
        self.type = type 
        self.pair = [labels.index(n) for n in pair]
        self.oth = [i for i in range(4) if i not in self.pair]
        self.loss = 0
    
    def unrooted_quartet(self):
        '''
        Returns the unrooted quartet, which is the same (self.pair | self.oth) regardless of the type.
        '''
        names = [self.labels[idx] for idx in self.pair] + [self.labels[idx] for idx in self.oth]
        if self.type == "BAL":
            return f"(({names[0]},{names[1]}),({names[2]},{names[3]}));"
        else:
            return f"({names[0]},({names[1]},({names[2]},{names[3]})));"

    def rooted_quartet(self, r = 'r'):
        '''
        Returns the rooted quartet in unrooted form by adding an extra node r on the root. Make sure r is not a name of any leaf.
        @param 
        r: the name of the extra node

        NOTE: When running ASTRAL with extra bipartitions, remember to add r to the bipartitions too!
        '''
        names = [self.labels[idx] for idx in self.pair] + [self.labels[idx] for idx in self.oth]
        if self.type == "BAL":
            return f"({r},({names[0]},{names[1]}),({names[2]},{names[3]}));"
        else:
            return f"({r},{names[0]},({names[1]},({names[2]},{names[3]})));"
    
    def get_unrooted_tuple(self):
        names = [self.labels[idx] for idx in self.pair] + [self.labels[idx] for idx in self.oth]
        return tuple(names)

    def add_states(self, l): 
        # l: list of names of leaves for which some state is present
        assert(all([n in self.labels for n in l]))
        if(len(l) == 0 or len(l) == 1 or len(l) == 4):
            # trivial cases: if no leaves, all leaves, or just one leaf has some state, then no losses are needed.
            return
        l = sorted([self.labels.index(n) for n in l])
        if(len(l) == 3):
            # BAL: born at root, loss at the one leaf without the state
            # LOP: no loss only if the states form a clade.
            if self.type == "LOP" and self.pair[0] not in l:
                # (a, (b, (c, d))) and l = [b, c, d]
                return 
            else:
                self.loss += 1
        else: 
            # len(l) = 2
            # no loss if states form a clade. 
            if self.type == "BAL":
                if sorted(self.pair) == l or sorted(self.oth) == l:
                    # ((a, b), (c, d)) and l = [a, b] or l = [c, d]
                    return 
                else:
                    self.loss += 2
            else:
                if sorted(self.oth) == l:
                    # no loss if clade
                    # (a, (b, (c, d))) and l = [c, d]
                    return 
                elif self.pair[1] in l:
                    # if (a, (b, (c, d))) and b is selected, only 1 loss is required
                    # [a, b], [b, c], [b, d]
                    self.loss += 1
                else:
                    # otherwise [a, c] or [a, d], and requires 2
                    self.loss += 2

def get_best_quartets_loss(data):
    '''
    Gets the best quartet trees from data
    @param: data is a dictionary name -> list of exhibited states. 
    @returns a list of RootedQuartets
    ''' 
    assert(len(data) == 4)
    names = list(data.keys())

    state_leaves = {}

    # get a dict state_leaves mapping state -> list of leaves which have that state
    for leaf, states in data.items():
        for s in states:
            if s not in state_leaves:
                state_leaves[s] = []
            state_leaves[s].append(leaf)
        
    # get all 15 possible trees
    possible_trees = [RootedQuartet("LOP", names, (names[i], names[j])) for i in range(4) for j in range(4) if i != j] + [RootedQuartet("BAL", names, (names[0], names[i])) for i in range(1, 4)]

    # calculate loss for all trees
    for _, leaves in state_leaves.items():
        for tree in possible_trees:
            tree.add_states(leaves)

    # rank and find trees with best score
    tree_scores = sorted([tree for tree in possible_trees], key=lambda x: x.loss) # sort in increasing order of loss
    best_trees = [t for t in tree_scores if t.loss == tree_scores[0].loss]

    return best_trees
