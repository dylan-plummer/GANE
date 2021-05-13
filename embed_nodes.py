import os
import numpy as np
import pandas as pd
import networkx as nx
import seaborn as sns
import matplotlib.pyplot as plt

from itertools import combinations
from scipy.sparse import coo_matrix
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from karateclub import NNSED, SymmNMF, DANMF
from karateclub import ASNE, FeatherNode, SINE, BANE, AE, MUSAE, TENE, TADW, FSCNMF

from cluster_core_utils import *

dataset = 'gavin'
out_dir = 'plots/%s' % dataset
os.makedirs(out_dir, exist_ok=True)

A = np.loadtxt('matrix/Attribute_%s.txt' % dataset)
M = np.loadtxt('matrix/Network_%s.txt' % dataset)

protein_df = pd.read_csv('dataset/%s_attr_vector.txt' % dataset, sep='\s+| |\t', header=None)
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

# get attribute matrix and convert to sparse format for attributed embedding methods
plt.matshow(nx.normalized_laplacian_matrix(G).A, cmap='prism')
plt.savefig(os.path.join(out_dir, 'attribute_matrix.png'))
plt.close()
attr_m = coo_matrix(A)
for method_name, method, community_detector in zip(['TENE', 'TADW', 'FSCNMF', 'ASNE', 'DANMF', 'SINE', 'BANE', 'AE', 'MUSAE'],  # the name of the embedding/community detection method
                                                   [TENE, TADW, FSCNMF, ASNE, DANMF, SINE, BANE, AE, MUSAE],  # the model class
                                                   [False, False, False, False, True, False, False, False, False]):  # whether the model can predict community membership
    print('Testing', method_name)
    os.makedirs(os.path.join(out_dir, method_name), exist_ok=True)
    try:
        model = method(workers=20)
    except TypeError:
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
    dist = pdist(z, metric='cosine')
    out_arr = np.vstack([comb, dist]).T
    out_df = pd.DataFrame(data=out_arr, columns=['p1', 'p2', 'edge', 'sim'])
    print(out_df)
    print('Removing non-edges')
    out_df = out_df.loc[out_df['edge'], ['p1', 'p2', 'sim']]
    print(out_df)
    out_df.to_csv(os.path.join(out_dir, method_name, '%s_sim.txt' % dataset), sep=' ', header=False, index=False)

    results_df = pd.DataFrame.copy(protein_df, deep=True)
    columns_names = ['protein'] + list(np.arange(0, len(results_df.columns) - 1))
    results_df.columns = [str(v) for v in columns_names]
    print(results_df.columns)
    z_tsne = TSNE(2).fit_transform(z)
    results_df['tsne_1'] = z_tsne[..., 0]
    results_df['tsne_2'] = z_tsne[..., 1]
    results_df['degree'] = degrees
    sns.scatterplot(data=results_df, x='tsne_1', y='tsne_2', hue='1', legend='brief', palette='Spectral')
    plt.title(method_name)
    plt.savefig(os.path.join(out_dir, method_name, '%s_tsne.png' % dataset))
    plt.close()
    print(results_df)
    results_df.to_csv(os.path.join(out_dir, method_name, '%s_embeddings.csv' % dataset))




