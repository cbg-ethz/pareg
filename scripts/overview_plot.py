import seaborn as sns
import matplotlib.pyplot as plt


def main(input_dirs, out_fname):
    plt.figure()

    sns.distplot([1,2,2,3,3,3,4])

    plt.savefig(out_fname)


if __name__ == '__main__':
    main(snakemake.input.input_dirs, snakemake.output.file)
