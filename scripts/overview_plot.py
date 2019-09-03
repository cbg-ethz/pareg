import os
import json
import datetime
import collections

import numpy as np
import pandas as pd

from dfply import *

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

        _, tool, group, source = dir_[:-1].split('/')
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


def annotate_correlation(x, y, *args, method='spearman', **kwargs):
    """Plot correlation.

    Adapted from https://github.com/mwaskom/seaborn/issues/1444
    """
    # compute correlation
    corr_r = x.corr(y, method)
    corr_text = f'{corr_r:2.2f}'.replace('0.', '.')

    # visualize correlation
    ax = plt.gca()
    ax.set_axis_off()

    if len(ax.collections) > 0:
        # TODO: handle this gracefully
        return

    marker_size = abs(corr_r) * 10000
    ax.scatter(
        .5, .5, marker_size, corr_r, alpha=0.6,
        cmap='vlag_r', vmin=-1, vmax=1,  # bwr_r
        transform=ax.transAxes)

    ax.annotate(
        corr_text,
        [.5, .5], xycoords='axes fraction',
        ha='center', va='center', fontsize=20)


def gentle_limit(limit_func, data, offset_scale=.05):
    """Set limits without cutting of data."""
    min_ = data.min()
    max_ = data.max()
    range_ = abs(max_ - min_)

    offset = range_ * offset_scale
    limit_func(min_ - offset, max_ + offset)


def format_axis(*args, ax=None, **kwargs):
    """Make axes look more minimalistic and fitting."""
    ax = ax or plt.gca()

    # sensible limits
    func_list = (ax.set_xlim, ax.set_ylim)
    assert len(func_list) >= len(args), 'woops'
    for func, data in zip(func_list, args):
        gentle_limit(func, data)

    # remove various annotations
    ax.tick_params(
        axis='both', which='both',
        left=False, labelleft=False,
        bottom=False, labelbottom=False)


def pvalue_scatterplots(df, out_dir):
    # prepare data
    df['idx'] = df['source'] + '|' + df['term']
    df_wide = (df.pivot(index='idx', columns='tool', values='trans_p_value')
                 .reset_index() >>
                 separate('idx', ['source', 'term'], sep=r'\|'))

    df_wide = df_wide.reindex(
        np.r_[
            df.loc[df['tool'].str.lower().argsort(), 'tool'].unique(),
            ['source', 'term']],
        axis=1)

    # plot
    with sns.plotting_context('talk', font_scale=2):
        g = sns.PairGrid(df_wide.dropna(), diag_sharey=False, height=5)

        g.map_upper(annotate_correlation, method='spearman')
        g.map_diag(
            sns.distplot, kde=False,
            bins=np.linspace(
                df_wide.drop(['source', 'term'], axis=1).values.ravel().min(),
                df_wide.drop(['source', 'term'], axis=1).values.ravel().max(),
                50))
        g.map_lower(sns.scatterplot, rasterized=True)

        g.map_lower(format_axis)
        g.map_diag(format_axis)  # [format_axis(ax=ax) for ax in g.diag_axes] # TODO: why does this not work?

        g.savefig(os.path.join(out_dir, 'pvalue_scatterplots.pdf'))


def runtime_overview(input_dirs, out_dir):
    # read data
    tmp = []
    for dir_ in input_dirs:
        fname = os.path.join(dir_, 'meta.json')
        with open(fname) as fd:
            meta_data = json.load(fd)

        _, tool, group, source = dir_[:-1].split('/')

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
    # organize input (TODO: do this with Snakemake)
    inp_data = collections.defaultdict(list)
    for dir_ in snakemake.input['input_dirs']:
        _, _, group, _, _ = dir_.split('/')
        inp_data[group].append(dir_)

    # execution
    for group, dir_list in inp_data.items():
        target_dir = os.path.join(snakemake.output['out_dir'], group)
        os.makedirs(target_dir, exist_ok=True)

        main(dir_list, target_dir)
