from utils import *
import pandas as pd
import numpy as np
import os,glob

#read in all file names to index
folder_path = 'data/'
len(glob.glob(os.path.join(folder_path, '*.tree')))

k = len(glob.glob(os.path.join(folder_path, '*.tree'))) #number of files in directory

c = np.empty(k, dtype=object) #empty numpy array for strings ------- check fixed lenght! faster dtype=char
index = np.empty(k, dtype=object)
i=0

for filename in glob.glob(os.path.join(folder_path, '*.tree')):
  with open(filename, 'r') as file:
    t = file.readlines()[6]
    index[i] = filename # get name to save string i.e. HumanBeta_Site17
    c[i] = t
    i= i+1


#main algorithm
w=5
k=len(index)
j=0
h1 = 0
h2 = 0
s = 0
data = np.empty([k, 5])
K = 2

for i in index:
    tree = Phylo.read(i, "nexus")
    # preprocessing
    collapse_tips(tree)
    add_names(tree)
    state_collapse(tree,K) 
    # to digraph 
    D = phylo_to_digraph(tree)
    x,e = low_freq_var(D,w)
    if x == 1:
        D.remove_nodes_from(e)
    c = conserved(D)
    if c == False:
        l = list(D.successors('0')) # find ancestors for each lineage
        p1 = all_states_paths(D,[l[0]]) # find paths by lineage
        p2 = all_states_paths(D,[l[1]]) # find paths by lineage
        p = [p1,p2]
        h1 = inter_homoplasy(D,p)
        h2 = intra_homoplasy(D,p)
        s = stepwise(D,p)
    data[j] = [x,c,h1,h2,s]
    j = j+1
        

  #results to dataframe and save to csv
df = pd.DataFrame(data=data, index=index)
df.to_csv(r'results.csv', header = False)
