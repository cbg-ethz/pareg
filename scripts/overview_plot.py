import os
import json

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


def read_result_data(input_dirs):
    df_list = []
    for dir_ in input_dirs:
        df = pd.read_csv(os.path.join(dir_, 'result.csv'))

        _, tool, source = dir_[:-1].split('/')
        df['tool'] = tool
        df['source'] = source

        df_list.append(df)

    return pd.concat(df_list)


def pvalue_overview(input_dirs, out_dir):
    df = read_result_data(input_dirs)

    g = sns.FacetGrid(df, col='source', row='tool')
    g.map(sns.distplot, 'pvalue', bins=100, kde=False)
    g.set(xlim=(0, 1))
    g.savefig(os.path.join(out_dir, 'overview.pdf'))


def runtime_overview(input_dirs, out_dir):
    # read data
    tmp = []
    for dir_ in input_dirs:
        fname = os.path.join(dir_, 'meta.json')
        with open(fname) as fd:
            meta_data = json.load(fd)

        _, tool, source = dir_[:-1].split('/')

        tmp.append({
            'tool': tool,
            'source': source,
            'runtime': meta_data['exec_time']
        })
    df = pd.DataFrame(tmp)

    # plot result
    plt.figure()

    sns.barplot(x='tool', y='runtime', hue='source', data=df)

    plt.ylabel('Runtime [s]')

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'runtime.pdf'))


def main(input_dirs, out_dir):
    pvalue_overview(input_dirs, out_dir)
    runtime_overview(input_dirs, out_dir)


if __name__ == '__main__':
    main(snakemake.input.input_dirs, snakemake.output.file)
