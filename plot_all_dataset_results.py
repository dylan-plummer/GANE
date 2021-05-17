import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn')

out_dir = 'plots'
os.makedirs(out_dir, exist_ok=True)

datasets = ['gavin', 'krogan2006core', 'dip', 'collins', 'biogrid']
df = pd.DataFrame()

for dataset in datasets:
    tmp_df = pd.read_csv('complex_results/complex_identification_results_%s.csv' % dataset)
    tmp_df['use_GO'] = tmp_df['use_GO'].apply(lambda v: True if v == 1 else False)
    tmp_df['dataset'] = dataset
    df = pd.concat([df, tmp_df])

df['composite'] = (df['Acc'] + df['F1']) * 0.5

methods = df['method'].unique()
method_map = {}
for i, method in enumerate(methods):
    method_map[method] = i
df['hue'] = df['method'].apply(lambda m: method_map[m])
print(df)
row = {'method': [], 'precision': [], 'recall': [], 'F1': [], 'Sn': [], 'PPV': [], 'Acc': [], 'composite': [], 'dataset': []}
for method in methods:
    for dataset in datasets:
        row['method'] += [method]
        row['dataset'] += [dataset]
        for metric in ['precision', 'recall', 'F1', 'Sn', 'PPV', 'Acc', 'composite']:
            try:
                tmp_df = df[(df['dataset'] == dataset) & (df['method'] == method)]
                go_val = tmp_df.iloc[0][metric]
                nogo_val = tmp_df.iloc[1][metric]
                row[metric] += [go_val - nogo_val]
            except IndexError:
                row[metric] += [0]
                pass

go_df = pd.DataFrame.from_dict(row)
print(go_df)
print(go_df['precision'])

for key in ['precision', 'recall', 'F1', 'Sn', 'PPV', 'Acc', 'composite']:
    g = sns.catplot(data=go_df, x='method', y=key, hue='dataset', kind='bar')
    #g.legend_.remove()
    plt.savefig(os.path.join(out_dir, '%s_go_change.png' % (key)))
    plt.close()

    g = sns.catplot(data=go_df, x='dataset', y=key,
                    hue='method', kind='bar', palette='tab10')
    # g.legend_.remove()
    plt.savefig(os.path.join(out_dir, '%s_dataset_go_change.png' % (key)))
    plt.close()

for key in ['precision', 'recall', 'F1', 'Sn', 'PPV', 'Acc', 'composite']:
    g = sns.catplot(data=df, x='method', y=key, order=df.drop_duplicates(subset='method').sort_values(by=key).method, hue='dataset', col='use_GO', kind='bar')
    #g.legend_.remove()
    plt.ylim(0, 1)
    plt.savefig(os.path.join(out_dir, '%s.png' % (key)))
    plt.close()

    g = sns.catplot(data=df, x='dataset', y=key,
                    hue='use_GO', col='method', kind='bar')
    # g.legend_.remove()
    plt.ylim(0, 1)
    plt.savefig(os.path.join(out_dir, '%s_method.png' % (key)))
    plt.close()

    g = sns.catplot(data=df, x='dataset', y=key,
                    hue='method', col='use_GO', kind='bar', palette='tab10')
    # g.legend_.remove()
    plt.ylim(0, 1)
    plt.savefig(os.path.join(out_dir, '%s_dataset.png' % (key)))
    plt.close()

