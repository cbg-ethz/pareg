import os
import json

import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


sns.set_context('talk')


def read_result_data(input_dirs):
    df_list = []
    for dir_ in input_dirs:
        df = pd.read_csv(os.path.join(dir_, 'result.csv'))

        _, tool, source = dir_[:-1].split('/')
        df['tool'] = tool
        df['source'] = source

        df_list.append(df)

    return pd.concat(df_list, ignore_index=True)


def pvalue_histograms(input_dirs, out_dir):
    df = read_result_data(input_dirs)

    g = sns.FacetGrid(df, col='source', row='tool', height=5)
    g.map(sns.distplot, 'p_value', bins=100, kde=False)
    g.set(xlim=(0, 1))
    g.savefig(os.path.join(out_dir, 'pvalue_histograms.pdf'))


def pvalue_scatterplots(input_dirs, out_dir):
    df = read_result_data(input_dirs)

    # custom pivot (TODO: make this better)
    df_piv = df.pivot_table(
        index='source', columns='tool', values='p_value', aggfunc=list)

    tmp = []
    for row in df_piv.itertuples():
        cur = []
        for col in df_piv.columns:
            cur.append([(col, v) for v in row._asdict()[col]])
        for foo in zip(*cur):
            t = {col: val for col, val in foo}
            t['source'] = row.Index
            tmp.append(t)
    df_wide = pd.DataFrame(tmp)

    # plot
    g = sns.PairGrid(df_wide.dropna(), hue='source', height=5)
    g = g.map_diag(sns.distplot, kde=False)
    g = g.map_offdiag(sns.scatterplot)
    g = g.add_legend()
    g.savefig(os.path.join(out_dir, 'pvalue_scatterplots.pdf'))


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
    pvalue_histograms(input_dirs, out_dir)
    pvalue_scatterplots(input_dirs, out_dir)
    runtime_overview(input_dirs, out_dir)


if __name__ == '__main__':
    main(snakemake.input.input_dirs, snakemake.output.file)
