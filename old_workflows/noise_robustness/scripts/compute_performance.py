from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.metrics import (
    precision_recall_curve, roc_curve, auc,
    average_precision_score, f1_score,
    confusion_matrix
)

import seaborn as sns
import matplotlib.pyplot as plt


def compute_pr(df, threshold=0.05):
    df['predicted_label'] = df['prediction_pvalue'] <= threshold

    tp = df[df['true_label'] & df['predicted_label']].shape[0]
    fp = df[~df['true_label'] & df['predicted_label']].shape[0]
    fn = df[df['true_label'] & ~df['predicted_label']].shape[0]

    try:
        precision = tp / (tp + fp)  # percentage of results which are relevant
    except ZeroDivisionError:
        precision = np.nan

    recall = tp / (tp + fn)  # percentage of all results which were classified

    return precision, recall


def main(fname_result, fname_termdatabase, outdir):
    outdir.mkdir()

    # load data
    df_result = pd.read_csv(fname_result)
    df_terms = pd.read_csv(fname_termdatabase)

    df = pd.DataFrame({
        'prediction_pvalue': df_result.set_index('term')['p_value'],
        'true_label': df_terms.drop_duplicates(subset=['term']).set_index('term')['term_is_enriched']
    })

    # compute measures
    df_pr = pd.DataFrame([
        (t, *compute_pr(df, t)) for t in np.logspace(-20, 0, 100)
    ], columns=['threshold', 'precision', 'recall'])

    # precision_list, recall_list, pr_thresholds = precision_recall_curve(df['true_label'], 1 - df['prediction_pvalue'])
    pr_auc = auc(df_pr['recall'], df_pr['precision'])

    # initialize plots
    s = 2
    fig_pr, ax_pr = plt.subplots(figsize=(s * 8, s * 6))

    # add to plot
    # ax_pr.plot(
    #     recall_list, precision_list,
    #     '-o',
    #     label=f'AUC: {pr_auc:.2}'
    # )
    sns.lineplot(
        data=df_pr, x='recall', y='precision',
        marker='o', estimator=None, sort=False,
        ax=ax_pr,
        label=f'AUC: {pr_auc:.2}'
    )

    # finalize plots
    ax_pr.set_xlabel('Recall')
    ax_pr.set_ylabel('Precision')
    ax_pr.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    fig_pr.tight_layout()
    fig_pr.savefig(outdir / f'pr_curve.pdf')


if __name__ == '__main__':
    main(
        snakemake.input.fname_result, snakemake.input.fname_termdatabase,
        Path(snakemake.output.outdir)
    )
