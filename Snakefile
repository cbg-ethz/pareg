import itertools


###
# setup

tool_list, = glob_wildcards('definitions/tools/{tool}/execute.py')

workdir: 'pipeline_run'

###
# rule definitions

def aggregate_input(wildcards):
    """Generate run-path for each tool/group/dataset."""
    checkpoint_output = checkpoints.provide_data.get(**wildcards).output.out_dir

    match = glob_wildcards(
        os.path.join(checkpoint_output, '{group}', '{dataset,[^./]+}/'))

    return [os.path.join('runs/', *p) + '/'
            for p in itertools.product(tool_list, [os.path.join(*p)
                                                   for p in zip(*match)])]


rule all:
    input:
        'overview_plots/',
        'performance/'

checkpoint provide_data:
    output:
        out_dir = directory('data/')
    shell:
        """
        # artifical data
        python3 {workflow.basedir}/scripts/generate_artificial_data.py {output.out_dir}

        # real data
        Rscript {workflow.basedir}/scripts/load_benchmark_data.R {output.out_dir}/real_data/
        """

rule execute:
    input:
        script = srcdir('definitions/tools/{tool}/execute.py'),
        data = 'data/{group}/{dataset}/'
    output:
        directory('runs/{tool}/{group}/{dataset}/')
    shell:
        """
        work_dir=$(pwd)
        inc_dir=$(cd {workflow.basedir}/scripts && pwd)
        cd {output}

        export PYTHONPATH="$inc_dir:$PYTHONPATH"
        python3 {input.script} \
            $work_dir/{input.data}/input.csv \
            $work_dir/{input.data}/terms.csv \
            $work_dir/{input.data}/term_network.csv
        """

rule summarize:
    input:
        input_dirs = aggregate_input
    output:
        out_dir = directory('overview_plots/')
    script:
        'scripts/overview_plot.py'

rule performance_evaluation:
    input:
        input_dirs = aggregate_input,
        data_dir = 'data/'
    output:
        out_dir = directory('performance/')
    script:
        'scripts/compute_performance.py'
