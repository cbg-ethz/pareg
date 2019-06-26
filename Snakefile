###
# setup

tool_list, = glob_wildcards('definitions/tools/{tool}/execute.py')

workdir: 'pipeline_run'

###
# rule definitions

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.provide_data.get(**wildcards).output.out_dir
    return expand(
        'runs/{tool}/{dataset}/',
        tool=tool_list,
        dataset=glob_wildcards(
            os.path.join(checkpoint_output, '{dataset,[^./]+}/')).dataset)

rule all:
    input:
        'overview_plots/'

checkpoint provide_data:
    output:
        out_dir = directory('data/')
    script:
        'scripts/load_benchmark_data.R'

rule execute:
    input:
        script = srcdir('definitions/tools/{tool}/execute.py'),
        data = 'data/{dataset}/'
    output:
        directory('runs/{tool}/{dataset}/')
    shell:
        """
        work_dir=$(pwd)
        inc_dir=$(cd {workflow.basedir}/scripts && pwd)
        cd {output}

        export PYTHONPATH="$inc_dir:$PYTHONPATH"
        python3 {input.script} \
            $work_dir/{input.data}/input.csv \
            $work_dir/{input.data}/terms.csv
        """

rule summarize:
    input:
        input_dirs = aggregate_input
    output:
        file = directory('overview_plots/')
    script:
        'scripts/overview_plot.py'
