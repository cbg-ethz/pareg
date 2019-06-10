import os
import sys
import glob
import shutil
import pandas as pd

import sh

from utils import Executor


class MyExecutor(Executor):
    def setup(self):
        # create ranked gene input
        self.gene_rank_file = 'genes.rnk'
        self.df_inp.sort_values('pvalue').to_csv(
            self.gene_rank_file, sep='\t', index=False, header=False)

        # store terms as gmt file (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29)
        self.term_file = 'terms.gmt'

        with open(self.term_file, 'w') as fd:
            for name, group in self.df_terms.groupby('name'):
                line = name + '\t'
                line += 'dummy_description\t'
                line += '\t'.join(group['gene'])

                fd.write(line + '\n')

    def execute(self):
        out_dir = 'output'
        shutil.rmtree(out_dir, ignore_errors=True)

        # run analysis
        root = os.path.dirname(os.path.realpath(__file__))
        sh.java(
            '-cp', os.path.join(root, 'gsea-3.0.jar'),
            '-Xmx512m', 'xtools.gsea.GseaPreranked',
            '-gmx', self.term_file,
            '-norm', None, '-nperm', 1000, '-scoring_scheme', 'classic',
            '-rpt_label', 'my_analysis',
            '-create_svgs', False, '-make_sets', True, '-plot_top_x', 20,
            '-rnd_seed', 42,
            '-set_max', 9999999, '-set_min', 0,
            '-zip_report', False,
            '-out', out_dir,
            '-gui', False,
            '-rnk', self.gene_rank_file)

        # parse output
        tmp = glob.glob(f'{out_dir}/*')
        assert len(tmp) == 1
        result_dir = tmp[0]

        idx = result_dir.split('.')[-1]
        fname_up = os.path.join(result_dir, f'gsea_report_for_na_pos_{idx}.xls')
        fname_down = os.path.join(result_dir, f'gsea_report_for_na_neg_{idx}.xls')

        df_up = pd.read_csv(fname_up, sep='\t')
        df_down = pd.read_csv(fname_down, sep='\t')

        # df_up['type'] = 'up'
        # df_down['type'] = 'down'

        self.df_result = pd.concat([
            df_up[['NAME', 'NOM p-val']],
            df_down[['NAME', 'NOM p-val']]
        ]).rename(columns={
            'NAME': 'term',
            'NOM p-val': 'p_value'
        })


if __name__ == '__main__':
    ex = MyExecutor(sys.argv[1], sys.argv[2])
    ex.run()
