import os
import math
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

plt.style.use('seaborn-dark')

dset_file = 'collins'
max_go = 32

for method in ['TENE', 'TADW', 'FSCNMF', 'ASNE', 'DANMF', 'SINE', 'BANE', 'AE', 'MUSAE']:
    for use_go in [True, False]:
        if use_go:
            dataset = dset_file + '_GO'
        else:
            dataset = dset_file
        try:
            df = pd.read_csv('plots/%s/%s/%s_embeddings.csv' % (dataset, method, dataset))
        except FileNotFoundError:
            continue
        print(df)
        #df['embedding'] = df['embedding'].apply(lambda s: [float(v) for v in s])
        n_go = len(df.columns) - 3
        if n_go > max_go:
            n_go = max_go
        nrows = int(math.sqrt(n_go))
        ncols = int(n_go / nrows)

        fig, axs = plt.subplots(nrows, ncols, figsize=(15, 15))

        for i, ax in enumerate(axs.reshape(-1)):
            print(i)
            if i == 0:
                sns.scatterplot(data=df, x='tsne_1', y='tsne_2', hue='degree', palette='viridis', alpha=0.75, ax=ax)
                ax.set_title('degree')
            else:
                sns.scatterplot(data=df, x='tsne_1', y='tsne_2', hue=str(i - 1), size='degree', palette='Spectral', alpha=0.8,
                                ax=ax, legend=False)
                ax.set_title(i)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel('')
            ax.set_ylabel('')

        #fig.subplots_adjust(wspace=0.01, hspace=0.01)
        fig.tight_layout()
        plt.savefig('plots/%s/%s.png' % (dataset, method))
        plt.close()