"""Adapted from GSEABenchmarkeR."""

import os

import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt


sns.set_context('talk')


def compute_relevance_score(df_obs, df_true):
    """Compute relevance score.

    1) Rank enrichment p-values (high rank -> top of GSEA ranking)
    2) Multiply ranks with corresponding relevance scores
    3) Sum result
    """
    # ranks = 1 - df_obs['p_value'].apply(
    #     lambda x: (df_obs['p_value'] <= x).mean())
    # return (ranks * df_true['relevance.score']).sum()

    tmp_obs = df_obs.set_index('term')
    tmp_true = df_true.set_index('term')

    ranks = 1 - tmp_obs['p_value'].apply(
        lambda x: (tmp_obs['p_value'] <= x).mean())
    return ranks.multiply(tmp_true['relevance.score']).sum()


def compute_random_score(df_true, perm=100):
    """Compute random relevance score distribution.

    Use this to assess significance of observed relevance score,
    i.e. how likely it is to observe a relevance score >= than the one obtained
    """
    # shuffled pvalue doesn't matter as only the rank is used
    return [compute_relevance_score(pd.DataFrame({
        'term': df_true['term'].sample(frac=1),
        'p_value': np.linspace(0, 1, df_true.shape[0])
    }), df_true) for _ in range(perm)]


def compute_optimal_score(df_true):
    """Compute best score one could obtain.

    Allows to compare score between multiple datasets
    (#geneset differs between phenotypes)

    Apply `compute_relevance_score` to case where GSEA ranking
    is identical to relevance score ranking
    """
    df_opt = pd.DataFrame({
        'term': df_true['term'],
        'p_value': np.linspace(0, 1, df_true.shape[0])
    })
    return compute_relevance_score(df_opt, df_true)


def compute_relative_relevance_score(df_obs, df_true):
    """Compute percentage of optimal score."""
    obs_score = compute_relevance_score(df_obs, df_true)
    opt_score = compute_optimal_score(df_true)

    return (obs_score / opt_score) * 100

def test():
    ea_ranks = pd.DataFrame({
        'term': ["a_bla", "d_hmm", "b_blu", "c_ha"],
        'p_value': [.01, .8, 0.2, 0.0005]
    })
    rel_ranks = pd.DataFrame({
        'term': ["a_bla", "b_blu", "c_ha", "d_hmm"],
        'relevance.score': [4, 3, 2, 1]
    })

    obs_score = compute_relevance_score(ea_ranks, rel_ranks)
    print(obs_score)

    rand_scores = compute_random_score(rel_ranks)
    print(rand_scores)

    res = (rand_scores >= obs_score).mean()
    print(res)

    opt_score = compute_optimal_score(rel_ranks)
    print(opt_score)

    opt_perc = compute_relative_relevance_score(ea_ranks, rel_ranks)
    print(opt_perc)

    import IPython; IPython.embed()


def main(input_dirs, data_dirs, out_dir):
    # compute scores
    tmp = []
    for dir_ in input_dirs:
        _, tool, dataset, _ = dir_.split('/')
        # read data
        fname_e = os.path.join(data_dirs, dataset, 'expected_terms.csv')
        df_expected = pd.read_csv(fname_e)

        fname_o = os.path.join(dir_, 'result.csv')
        df_observed = pd.read_csv(fname_o)

        # compute and store
        score = compute_relative_relevance_score(df_observed, df_expected)
        tmp.append({
            'tool': tool,
            'dataset': dataset,
            'score': score
        })

    df = pd.DataFrame(tmp)

    # plot result
    plt.figure(figsize=(8, 6))

    sns.boxplot(x='tool', y='score', data=df)
    sns.stripplot(x='tool', y='score', data=df)

    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'score_boxplot.pdf'))


if __name__ == '__main__':
    main(
        snakemake.input['input_dirs'], snakemake.input['data_dir'],
        snakemake.output['out_dir'])
