import os
import shutil

import numpy as np
import pandas as pd

from tqdm import trange


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
        'term': ['foo'] * 100,
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


def main(output_dir):
    sliding_window_data(output_dir)


if __name__ == '__main__':
    main('definitions/data/')
