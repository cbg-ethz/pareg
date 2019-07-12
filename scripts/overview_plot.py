import os
import json
import datetime

import numpy as np
import pandas as pd

from statsmodels.sandbox.stats.multicomp import multipletests

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

from tool_handlers import TRANSFORMER_DICT


sns.set_context('talk')


@FuncFormatter
def format_seconds(x, pos):
    return str(datetime.timedelta(seconds=round(x)))


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
    g.map(sns.distplot, 'p_value', kde=False)
    g.savefig(os.path.join(out_dir, 'pvalue_histograms.pdf'))


def annotate_correlation(*args, method, **kwargs):
    """Plot correlation.

    Adapted from https://github.com/mwaskom/seaborn/issues/1444
    """
    # compute correlation
    corr_r = args[0].corr(args[1], method)
    corr_text = f'{corr_r:2.2f}'.replace('0.', '.')

    # visualize correlation
    ax = plt.gca()
    ax.set_axis_off()

    marker_size = abs(corr_r) * 10000
    ax.scatter(
        .5, .5, marker_size, corr_r, alpha=0.6,
        cmap='vlag', vmin=-1, vmax=1,  # bwr_r
        transform=ax.transAxes)

    ax.annotate(
        corr_text,
        [.5, .5], xycoords='axes fraction',
        ha='center', va='center', fontsize=20)


def pvalue_scatterplots(df, out_dir):
    # custom pivot (TODO: make this better)
    df_piv = df.pivot_table(
        index='source', columns='tool', values='trans_p_value', aggfunc=list)

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
        np.r_[
            df.loc[df['tool'].str.lower().argsort(), 'tool'].unique(),
            ['source']],
        axis=1)

    # plot
    g = sns.PairGrid(df_wide.dropna(), height=5)  # , hue='source'
    g = g.map_upper(annotate_correlation, method='spearman')
    g = g.map_diag(sns.distplot, kde=False)
    g = g.map_lower(sns.scatterplot, rasterized=True)
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
    plt.ylabel('Runtime [h:m:s]')

    plt.xticks(rotation=90)
    plt.yscale('log')

    plt.gca().yaxis.set_major_formatter(format_seconds)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'runtime.pdf'))


def significant_term_counts(df, out_dir):
    # gather
    # df.groupby(['tool', 'source'])['p_value'].apply(lambda x: multipletests(x, method='fdr_bh')[1])

    df_sig = (
        df.groupby(['tool', 'source'])['p_value']
          .apply(
            lambda x: (TRANSFORMER_DICT[x.name[0]].threshold_pvalues(x).sum()))
          .reset_index()
          .rename(columns={'p_value': 'term_count'})
    )

    # assume that we only use one term database (for now...)
    assert df.groupby(['tool', 'source'])['term'].count().unique().size == 1
    total_term_count = df.groupby(['tool', 'source'])['term'].count().iloc[0]

    df_long = pd.melt(df_sig, id_vars=['tool', 'source'])
    df_long['value'] /= total_term_count

    # plot
    plt.figure(figsize=(8, 8))

    sns.boxplot(
        x='tool', y='value', data=df_long,
        order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique())
    sns.stripplot(
        x='tool', y='value', data=df_long,
        order=df.loc[df['tool'].str.lower().argsort(), 'tool'].unique())

    plt.xlabel('Tool')
    plt.ylabel('%terms with $\mathrm{pvalue} < 0.05$')

    plt.xticks(rotation=90)
    plt.ylim((-.1, 1.1))

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'term_counts.pdf'))


def main(input_dirs, out_dir):
    # read data
    df = read_result_data(input_dirs)
    df.to_csv(os.path.join(out_dir, 'results.csv'), index=False)

    # transform data
    df['trans_p_value'] = df.apply(
        lambda x: TRANSFORMER_DICT[x.tool].transform_pvalues([x.p_value])[0],
        axis=1)
    df['trans_p_value'] = (df.groupby(['tool', 'source'])['trans_p_value']
                             .transform(lambda x: x.fillna(x.max())))

    # plot data
    pvalue_histograms(df, out_dir)
    pvalue_scatterplots(df, out_dir)
    significant_term_counts(df, out_dir)

    runtime_overview(input_dirs, out_dir)


if __name__ == '__main__':
    main(snakemake.input.input_dirs, snakemake.output.file)
