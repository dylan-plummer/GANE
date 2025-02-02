import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('seaborn')

out_dir = 'plots'
dataset = 'gavin'
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv('complex_results/complex_identification_results_%s.csv' % dataset)
methods = df['method'].unique()
method_map = {}
for i, method in enumerate(methods):
    method_map[method] = i
df['hue'] = df['method'].apply(lambda m: method_map[m])
print(df)

for key in ['precision', 'recall', 'F1', 'Sn', 'PPV', 'Acc']:
    g = sns.barplot(data=df, x='method', y=key, order=df.drop_duplicates(subset='method').sort_values(by=key).method, hue='use_GO')
    #g.legend_.remove()
    plt.ylim(0, 1)
    plt.savefig(os.path.join(out_dir, '%s_%s.png' % (dataset, key)))
    plt.close()
