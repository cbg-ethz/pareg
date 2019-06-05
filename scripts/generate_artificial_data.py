import os
import shutil

import numpy as np
import pandas as pd

from tqdm import trange


def generate_input(i):
    window_size = 10

    genes = [f'g{i:02}' for i in range(90, 110)]
    p_values = np.concatenate([
        np.random.beta(1, 1, i),  # uniform
        np.random.beta(0.1, 1, window_size),  # "significant" p-values
        np.random.beta(1, 1, 20-i-window_size)  # uniform
    ])

    df = pd.DataFrame({
        'gene': genes,
        'pvalue': p_values
    })
    return df


def generate_terms():
    df = pd.DataFrame({
        'name': ['foo'] * 100,
        'gene': [f'g{i:02}' for i in range(100, 200)],
    })
    return df


def main(output_dir):
    for i in trange(10):
        # generate data
        df_genes = generate_input(i)
        df_terms = generate_terms()

        # store data
        target_dir = os.path.join(output_dir, f'artificial_data_{i:02}')
        shutil.rmtree(target_dir, ignore_errors=True)
        os.makedirs(target_dir)

        df_genes.to_csv(os.path.join(target_dir, 'input.csv'), index=False)
        df_terms.to_csv(os.path.join(target_dir, 'terms.csv'), index=False)


if __name__ == '__main__':
    main('definitions/data/')
