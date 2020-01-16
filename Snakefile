###
# setup

tool_list, = glob_wildcards('definitions/tools/{tool}/execute.py')
data_sources, = glob_wildcards('definitions/data_generators/{source}.R')

workdir: 'pipeline_run'


###
# rule definitions

rule all:
    input:
        expand('comparison/{source}/', source=data_sources)


rule generate_data:
    input:
        script = srcdir('definitions/data_generators/{source}.R')
    output:
        cell_fname = 'data/{source}/cell_matrix.csv',
        expr_fname = 'data/{source}/expression_matrix.csv',
        info_fname = 'data/{source}/condition_info.csv'
    shell:
        """
        Rscript "{input.script}" \
            "{output.cell_fname}" \
            "{output.expr_fname}" \
            "{output.info_fname}"
        """


rule perform_dea:
    input:
        expr_fname = 'data/{source}/expression_matrix.csv',
        info_fname = 'data/{source}/condition_info.csv'
    output:
        dea_fname = 'data/{source}/dea_result.csv'
    script:
        'scripts/perform_dea.R'


rule execute_tool:
    input:
        script = srcdir('definitions/tools/{tool}/execute.py'),
        cell_fname = 'data/{source}/cell_matrix.csv',
        expr_fname = 'data/{source}/expression_matrix.csv',
        info_fname = 'data/{source}/condition_info.csv',
        dea_fname = 'data/{source}/dea_result.csv'
    output:
        run_dir = directory('runs/{tool}/{source}/run_dir'),
        result_fname = 'runs/{tool}/{source}/enrichment_result.csv',
        meta_fname = 'runs/{tool}/{source}/meta.json'
    shell:
        """
        work_dir=$(pwd)
        inc_dir=$(cd {workflow.basedir}/scripts && pwd)

        export PYTHONPATH="$inc_dir:$PYTHONPATH"
        python3 {input.script} \
            "$work_dir/{output.run_dir}" \
            "$work_dir/{input.cell_fname}" \
            "$work_dir/{input.expr_fname}" \
            "$work_dir/{input.info_fname}" \
            "$work_dir/{input.dea_fname}" \
            "$work_dir/{output.result_fname}" \
            "$work_dir/{output.meta_fname}"
        """


rule compare_tools:
    input:
        expand('runs/{tool}/{{source}}/enrichment_result.csv', tool=tool_list)
    output:
        out_dir = directory('comparison/{source}/')
    log:
        notebook = 'notebooks/ToolComparison.{source}.ipynb'
    notebook:
        'notebooks/ToolComparison.ipynb'
