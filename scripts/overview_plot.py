import os
import json

import pandas as pd

from statsmodels.sandbox.stats.multicomp import multipletests

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


def pvalue_histograms(df, out_dir):
    g = sns.FacetGrid(
        df, col='source', row='tool',
        row_order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique(),
        sharex=False, sharey=False,
        height=5)
    g.map(sns.distplot, 'p_value', bins=100, kde=False)
    g.savefig(os.path.join(out_dir, 'pvalue_histograms.pdf'))


def pvalue_scatterplots(df, out_dir):
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
    df_wide = df_wide.reindex(
        df.loc[df['tool'].str.lower().argsort(), 'tool'].unique(), axis=1)

    # plot
    g = sns.PairGrid(df_wide.dropna(), height=5)  # , hue='source'
    g = g.map_diag(sns.distplot, kde=False)
    g = g.map_offdiag(sns.scatterplot, rasterized=True)
    # g = g.add_legend()
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
    plt.figure(figsize=(8, 8))

    sns.boxplot(
        x='tool', y='runtime', data=df,
        order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique())
    sns.stripplot(
        x='tool', y='runtime', data=df,
        order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique())

    plt.xlabel('Tool')
    plt.ylabel('Runtime [s]')

    plt.xticks(rotation=90)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'runtime.pdf'))


def significant_term_counts(df, out_dir):
    # gather
    df['p_value_corrected'] = multipletests(df['p_value'], method='fdr_bh')[1]

    df_sig = (df.groupby(['tool', 'source'])
                .aggregate(lambda x: (x <= .05).sum())
                .rename(columns={
                    'p_value': 'uncorrected',
                    'p_value_corrected': 'corrected'})
                .reset_index())

    df_long = pd.melt(df_sig, id_vars=['tool', 'source'])

    # plot
    plt.figure(figsize=(12, 8))

    sns.boxplot(
        x='tool', y='value', hue='variable', data=df_long,
        order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique())
    sns.stripplot(
        x='tool', y='value', hue='variable', data=df_long,
        order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique())

    plt.xlabel('Tool')
    plt.ylabel('|Terms with $\mathrm{pvalue} < 0.05$|')

    plt.xticks(rotation=90)

    handles, labels = plt.gca().get_legend_handles_labels()
    plt.legend(
        handles[:2], labels[:2],
        bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'term_counts.pdf'))


def main(input_dirs, out_dir):
    df = read_result_data(input_dirs)
    df.to_csv(os.path.join(out_dir, 'results.csv'), index=False)

    pvalue_histograms(df, out_dir)
    pvalue_scatterplots(df, out_dir)
    significant_term_counts(df, out_dir)

    runtime_overview(input_dirs, out_dir)


if __name__ == '__main__':
    main(snakemake.input.input_dirs, snakemake.output.file)
