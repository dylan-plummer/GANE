import os
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import combinations
from scipy.sparse import coo_matrix
from scipy.spatial.distance import pdist
from sklearn.manifold import TSNE
from karateclub import DANMF
from karateclub import ASNE, SINE, BANE, AE, MUSAE, TENE, TADW, FSCNMF

from cluster_core_utils import *

#dataset = 'biogrid'
dset_file = 'gavin'
dataset = dset_file
n_exp = 1
use_go = True
if use_go:
    dataset += '_GO'
out_dir = 'plots/%s' % dataset
os.makedirs(out_dir, exist_ok=True)

A = np.loadtxt('matrix/Attribute_%s.txt' % dset_file)
M = np.loadtxt('matrix/Network_%s.txt' % dset_file)

protein_df = pd.read_csv('dataset/%s_attr_vector.txt' % (dset_file), sep='\s+| |\t', header=None)
print(protein_df)
protein_list = np.array(protein_df.iloc[:, 0].values)
print(protein_list)
print(len(protein_list))
R,C = np.triu_indices(len(protein_list),1)
edge_mask = []
for i, j in combinations(np.arange(M.shape[0]), 2):
    edge_mask.append(M[i][j] == 1)
edge_mask = np.array(edge_mask)  # 1d boolean mask indicating presence of edges in flattened upper triangle
comb = np.vstack([protein_list[R], protein_list[C], edge_mask])  # this plus a column for similarities will be the final output

G = nx.convert_matrix.from_numpy_array(M)  # create graph from adjacency matrix
print(len(G.edges()))

degrees = []
for (_, d) in G.degree():
    degrees.append(d)

# plot degree distribution
degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
dmax = max(degree_sequence)
plt.loglog(degree_sequence, "b-", marker="o")
plt.title("Degree rank plot")
plt.ylabel("degree")
plt.xlabel("rank")
plt.savefig(os.path.join(out_dir, 'degree_distribution.png'))
plt.close()

plt.matshow(nx.normalized_laplacian_matrix(G).A, cmap='nipy_spectral')
plt.savefig(os.path.join(out_dir, 'laplacian_norm.png'))
plt.close()
# get attribute matrix and convert to sparse format for attributed embedding methods
if use_go:
    attr_m = coo_matrix(A)
else:
    attr_m = coo_matrix(np.zeros_like(A))

# dict for storing results (will be converted to DataFrame)
metrics = {'method': [], 'precision': [], 'recall': [], 'F1': [], 'Sn': [], 'PPV': [], 'Acc': []}
for method_name, method, community_detector in zip(['TENE', 'TADW', 'FSCNMF', 'ASNE', 'DANMF', 'SINE', 'BANE', 'AE', 'MUSAE'],  # the name of the embedding/community detection method
                                                   [TENE, TADW, FSCNMF, ASNE, DANMF, SINE, BANE, AE, MUSAE],  # the model class
                                                   [False, False, False, False, True, False, False, False, False]):  # whether the model can predict community membership
    for exp_i in range(n_exp):
        print('Testing', method_name, exp_i)
        os.makedirs(os.path.join(out_dir, method_name), exist_ok=True)
        model = method()
        model.allow_disjoint = True
        community_df = pd.DataFrame()
        try:
            if community_detector:  # model can predict community membership
                model.fit(G)
                z = model.get_embedding()
            else:  # otherwise run density-based clustering in embedding space
                model.fit(G, attr_m)
                z = model.get_embedding()
        except AssertionError as e:
            print(e)
            continue
        dist = pdist(z, metric='cosine')  # compute pairwise similarity between every pair of nodes
        out_arr = np.vstack([comb, dist]).T
        out_df = pd.DataFrame(data=out_arr, columns=['p1', 'p2', 'edge', 'sim'])
        out_df = out_df.loc[out_df['edge'], ['p1', 'p2', 'sim']]  # filter out rows where proteins did not share an edge
        out_df.to_csv(os.path.join(out_dir, method_name, '%s_sim.txt' % dataset), sep=' ', header=False, index=False)

        if exp_i == 0:
            results_df = pd.DataFrame.copy(protein_df, deep=True)
            columns_names = ['protein'] + list(np.arange(0, len(results_df.columns) - 1))
            results_df.columns = [str(v) for v in columns_names]
            z_tsne = TSNE(2).fit_transform(z)
            results_df['tsne_1'] = z_tsne[..., 0]
            results_df['tsne_2'] = z_tsne[..., 1]
            results_df['degree'] = degrees
            results_df.to_csv(os.path.join(out_dir, method_name, '%s_embeddings.csv' % dataset))

        cluster(dataset, dset_file, method_name)
        precision, recall, F1, Sn, PPV, Acc = eval(dataset)
        metrics['method'] += [method_name]
        metrics['precision'] += [precision]
        metrics['recall'] += [recall]
        metrics['F1'] += [F1]
        metrics['Sn'] += [Sn]
        metrics['PPV'] += [PPV]
        metrics['Acc'] += [Acc]

df = pd.DataFrame.from_dict(metrics)
df.to_csv('complex_identification_results.csv')

methods = df['method'].unique()
method_map = {}
for i, method in enumerate(methods):
    method_map[method] = i
df['hue'] = df['method'].apply(lambda m: method_map[m])
print(df)

for key in ['precision', 'recall', 'F1', 'Sn', 'PPV', 'Acc']:
    g = sns.barplot(data=df, x='method', y=key, order=df.sort_values(by=key).method)
    g.legend_.remove()
    plt.ylim(0, 1)
    plt.savefig(os.path.join(out_dir, '%s.png' % (key)))
    plt.close()




