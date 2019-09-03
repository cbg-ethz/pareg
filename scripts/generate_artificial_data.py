import os
import sys
import shutil

import numpy as np
import pandas as pd

from tqdm import tqdm, trange


def sliding_window_data(output_dir):
    """Generate input data whose significant genes slide over terms.

    1) uniform inside, beta outside
    2) half/half
    3) uniform outside, beta inside
    """
    # helper functions
    def generate_input(i):
        window_size = 10

        genes = [f'g{i:02}' for i in range(90, 110)]
        p_values = np.concatenate([
            np.random.beta(1, 1, size=i),  # uniform
            np.random.beta(0.1, 1, size=window_size),  # "significant" p-values
            np.random.beta(1, 1, size=20-i-window_size)  # uniform
        ])

        return pd.DataFrame({
            'gene': genes,
            'p_value': p_values
        })

    # generate data
    df_terms = pd.DataFrame({
        'term': 'foo',
        'gene': [f'g{i:02}' for i in range(100, 200)],
    })
    df_expected = pd.DataFrame({
        'term': ['foo'],
        'relevance.score': [10],
    })

    for i in trange(10):
        # generate data
        df_genes = generate_input(i)

        # store data
        target_dir = os.path.join(output_dir, f'artificial_data_{i:02}')
        shutil.rmtree(target_dir, ignore_errors=True)
        os.makedirs(target_dir)

        df_genes.to_csv(os.path.join(target_dir, 'input.csv'), index=False)
        df_terms.to_csv(os.path.join(target_dir, 'terms.csv'), index=False)
        df_expected.to_csv(os.path.join(target_dir, 'expected_terms.csv'), index=False)


def sanity_check(output_dir):
    """Generate artificial data.

    Dataset which shows how tools behave for more/fewer/equal DE gene counts.
    """
    # generate gene set
    genes_sig = [f'g{i:03}' for i in range(0, 300)]
    genes_uni = [f'g{i:03}' for i in range(300, 600)]
    genes_out = [f'g{i:03}' for i in range(600, 900)]

    df_genes = pd.DataFrame({
        'gene': np.r_[genes_sig, genes_uni],
        'p_value': np.r_[
            np.random.beta(0.1, 1, size=len(genes_sig)),
            np.random.beta(1, 1, size=len(genes_uni))
        ]
    })

    # generate terms
    N = 100
    frac_list = [
        [.5, .5, 0],  # sig, uni, out
        [.8, .2, 0],
        [.2, .8, 0]
    ]

    for i, fracs in enumerate(tqdm(frac_list)):
        assert sum(fracs) == 1

        gs_set = [genes_sig, genes_uni, genes_out]
        assert len(fracs) == len(gs_set)

        df_list = []
        for j in range(10):
            term_genes = np.concatenate([np.random.choice(
                gs, size=round(N*f), replace=False)
                                         for f, gs in zip(fracs, gs_set)])
            assert len(term_genes) == N

            tmp = pd.DataFrame({
                'term': f'gs_{j:02}',
                'gene': term_genes
            })
            df_list.append(tmp)
        df_terms = pd.concat(df_list)

        # dummy ranking expectations
        df_expected = pd.DataFrame({
            'term': ['gs_09', 'gs_01', 'gs_00'],
            'relevance.score': [10, 5, 2],
        })

        # save results
        target_dir = os.path.join(output_dir, f'sanity_check_{i:02}')
        shutil.rmtree(target_dir, ignore_errors=True)
        os.makedirs(target_dir)

        df_genes.to_csv(os.path.join(target_dir, 'input.csv'), index=False)
        df_terms.to_csv(os.path.join(target_dir, 'terms.csv'), index=False)
        df_expected.to_csv(os.path.join(target_dir, 'expected_terms.csv'), index=False)


def main(output_dir):
    data_generators = [sanity_check, scoring_check]
    for func in data_generators:
        name = func.__name__
        print(f' >> {name} <<')

        outdir = os.path.join(output_dir, name)
        func(outdir)


if __name__ == '__main__':
    main(sys.argv[1])
