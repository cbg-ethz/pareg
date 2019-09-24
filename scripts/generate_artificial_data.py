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


class DataGenerator:
    """Given input genes, create terms."""

    def __init__(self, input_size=1000):
        np.random.seed(42)
        self._setup(input_size)

    def _setup(self, input_size):
        """Initialize data.

        Create equally sized gene sets
        """
        step = input_size // 3
        self.genes_sig = [f'g{i:03}' for i in range(0, step)]
        self.genes_uni = [f'g{i:03}' for i in range(step, 2*step)]
        self.genes_out = [f'g{i:03}' for i in range(2*step, 3*step)]

        self.df_genes = pd.DataFrame({
            'gene': np.r_[self.genes_sig, self.genes_uni],
            'p_value': np.r_[
                np.random.beta(0.1, 1, size=len(self.genes_sig)),
                np.random.beta(1, 1, size=len(self.genes_uni))
            ]
        })

        print(f'Input gene set has size {self.df_genes.shape[0]}')

    def generate_term(self, frac_sig, frac_uni, frac_out, name, term_size=100):
        """Generate custom gene setself.

        The created gene set has a certain fraction of genes from input
        which have a significant and non-significant p-value, as well as
        genes which are in the input set.
        """
        # setup
        fracs = [frac_sig, frac_uni, frac_out]
        gs_set = [self.genes_sig, self.genes_uni, self.genes_out]

        # generate term
        term_genes = np.concatenate([np.random.choice(
                                        gs,
                                        size=np.round(term_size*f).astype(int),
                                        replace=False)
                                     for f, gs in zip(fracs, gs_set)])

        df_term = pd.DataFrame({
            'term': name,
            'gene': term_genes
        })

        # sanity checks
        assert sum(fracs) == 1, fracs
        assert len(fracs) == len(gs_set)
        assert df_term.shape[0] == term_size

        return df_term

    @staticmethod
    def generate_term_network(df_terms):
        """Generate adjacency matrix for term network.

        Given a concatenated term dataframe, generate a network with weights.
        """
        term_num = df_terms['term'].nunique()
        adja_mat = pd.DataFrame(np.ones(shape=(term_num, term_num)))
        np.fill_diagonal(adja_mat.values, 0)

        unique_terms = df_terms['term'].unique()
        adja_mat.index = unique_terms
        adja_mat.columns = unique_terms

        return adja_mat


def sanity_check(output_dir):
    """Generate artificial data.

    Dataset which shows how tools behave for more/fewer/equal DE gene counts.
    """
    # generate terms
    generator = DataGenerator(900)

    frac_list = [
        [.5, .5, 0],  # sig, uni, out
        [.8, .2, 0],
        [.2, .8, 0]
    ]

    for i, fracs in enumerate(tqdm(frac_list)):
        df_list = []
        for j in range(10):
            tmp = generator.generate_term(*fracs, f'gs_{j:02}', 100)
            df_list.append(tmp)
        df_terms = pd.concat(df_list)

        # term network
        term_network = generator.generate_term_network(df_terms)

        # dummy ranking expectations
        df_expected = pd.DataFrame({
            'term': ['gs_09', 'gs_01', 'gs_00'],
            'relevance.score': [10, 5, 2],
        })

        # save results
        target_dir = os.path.join(output_dir, f'sanity_check_{i:02}')
        shutil.rmtree(target_dir, ignore_errors=True)
        os.makedirs(target_dir)

        generator.df_genes.to_csv(
            os.path.join(target_dir, 'input.csv'), index=False)
        df_terms.to_csv(
            os.path.join(target_dir, 'terms.csv'), index=False)

        term_network.to_csv(
            os.path.join(target_dir, 'term_network.csv'))

        df_expected.to_csv(
            os.path.join(target_dir, 'expected_terms.csv'), index=False)


def scoring_check(output_dir):
    """Check final tool scoring scheme.

    Data with meaningful output rankings.
    """
    # generate terms
    generator = DataGenerator(900)

    frac_list = [[x, 1-x, 0]  # sig, uni, out
                 for x in np.arange(0.5, 1, .01)]

    for i in trange(10):
        df_list = []
        for j, fracs in enumerate(frac_list):
            tmp = generator.generate_term(*fracs, f'gs_{j:02}', 100)
            df_list.append(tmp)
        df_terms = pd.concat(df_list)

        # term network
        term_network = generator.generate_term_network(df_terms)

        # dummy ranking expectations
        df_expected = pd.DataFrame({
            'term': df_terms['term'].unique()[::-1],
            'relevance.score': [(i+1)*10
                                for i in reversed(range(df_terms['term'].unique().size))],
        })

        # save results
        target_dir = os.path.join(output_dir, f'scoring_check_{i:02}')
        shutil.rmtree(target_dir, ignore_errors=True)
        os.makedirs(target_dir)

        generator.df_genes.to_csv(
            os.path.join(target_dir, 'input.csv'), index=False)
        df_terms.to_csv(
            os.path.join(target_dir, 'terms.csv'), index=False)

        term_network.to_csv(
            os.path.join(target_dir, 'term_network.csv'))

        df_expected.to_csv(
            os.path.join(target_dir, 'expected_terms.csv'), index=False)


def robustness_check(output_dir):
    """Check robustness of tools for decreasing signal.

    Generate "(dis)similar" terms and vary input to reduce overlap.
    """
    # prepare data generation
    generator = DataGenerator(1000)

    for i, x in enumerate(tqdm(np.linspace(0.5, 1., 4))):
        for j in range(3):
            df_list = []

            # "master" terms
            for k in range(2):
                df_master = generator.generate_term(
                    frac_sig=.9, frac_uni=.1, frac_out=0,
                    name=f'master_{k:02}', term_size=100)
                df_list.append(df_master)

            # generate other terms
            df_connected = generator.generate_term(
                frac_sig=x, frac_uni=1-x, frac_out=0,
                name=f'connected', term_size=50)
            # df_isolated = generator.generate_term(
            #     frac_sig=x, frac_uni=1-x, frac_out=0,
            #     name=f'isolated', term_size=50)
            df_list.extend([df_connected])  # df_isolated

            df_terms = pd.concat(df_list)

            # term network
            term_network = generator.generate_term_network(df_terms)

            # dummy ranking expectations
            df_expected = pd.DataFrame({
                'term': df_terms['term'].unique()[::-1],
                'relevance.score': [(i+1)*10
                                    for i in reversed(range(df_terms['term'].unique().size))],
            })

            # save results
            target_dir = os.path.join(output_dir, f'robustness_check_{i:02}__{j:02}')
            shutil.rmtree(target_dir, ignore_errors=True)
            os.makedirs(target_dir)

            generator.df_genes.to_csv(
                os.path.join(target_dir, 'input.csv'), index=False)
            df_terms.to_csv(
                os.path.join(target_dir, 'terms.csv'), index=False)

            term_network.to_csv(
                os.path.join(target_dir, 'term_network.csv'))

            df_expected.to_csv(
                os.path.join(target_dir, 'expected_terms.csv'), index=False)


def main(output_dir):
    data_generators = [
        sanity_check,
        scoring_check,
        robustness_check
    ]

    for func in data_generators:
        name = func.__name__
        print(f' >> {name} <<')

        outdir = os.path.join(output_dir, name)
        func(outdir)


if __name__ == '__main__':
    main(sys.argv[1])
