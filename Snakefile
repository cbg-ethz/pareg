###
# setup

tool_list, = glob_wildcards('definitions/tools/{tool}/execute.py')
source_list, = glob_wildcards('definitions/data/{source}/terms.csv')

print(tool_list)
print(source_list)


###
# rule definitions

rule all:
    input:
        'results/overview.pdf'

rule execute:
    input:
        script = srcdir('definitions/tools/{tool}/execute.py'),
        data = srcdir('definitions/data/{source}/')
    output:
        directory('results/{tool}/{source}/')
    shell:
        """
        inc_dir=$(cd scripts && pwd)
        cd {output}

        export PYTHONPATH="$inc_dir:$PYTHONPATH"
        python3 {input.script} {input.data}/input.csv {input.data}/terms.csv
        """

rule summarize:
    input:
        input_dirs = expand(
            'results/{tool}/{source}/',
            tool=tool_list, source=source_list)
    output:
        file = 'results/overview.pdf'
    script:
        'scripts/overview_plot.py'
