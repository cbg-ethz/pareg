###
# setup

tool_list, = glob_wildcards('definitions/tools/{tool}/execute.sh')
source_list, = glob_wildcards('definitions/data/{source}/terms.csv')

print(source_list)

###
# rule definitions

rule all:
    input:
        'results/overview.pdf'

rule execute:
    input:
        script = srcdir('definitions/tools/{tool}/execute.sh'),
        data = 'definitions/data/{source}/'
    output:
        directory('results/{tool}/{source}/')
    shell:
        """
        cd {output}
        {input.script} {input.data}/input.csv {input.data}/terms.csv
        """

rule summarize:
    input:
        input_dirs = expand('results/{tool}/{source}/', tool=tool_list, source=source_list)
    output:
        file = 'results/overview.pdf'
    script:
        'scripts/overview_plot.py'
