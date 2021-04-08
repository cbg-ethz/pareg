###
# setup

tool_list, = glob_wildcards('definitions/tools/{tool}/execute.py')
data_sources, = glob_wildcards('definitions/data_generators/{source}.R')

workdir: 'pipeline_run'


###
# rule definitions

rule all:
    input:
        expand('comparison/{source}/', source=data_sources),
        expand('data/{source}/data_stats/', source=data_sources),
        expand('data/{source}/dea_stats/', source=data_sources)


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


rule analyze_expression_data:
    input:
        expr_fname = 'data/{source}/expression_matrix.csv',
        info_fname = 'data/{source}/condition_info.csv'
    output:
        out_dir = directory('data/{source}/data_stats/')
    log:
        notebook = 'notebooks/AnalyzeExpressionData.{source}.ipynb'
    notebook:
        'notebooks/AnalyzeExpressionData.ipynb'


rule perform_dea:
    input:
        expr_fname = 'data/{source}/expression_matrix.csv',
        info_fname = 'data/{source}/condition_info.csv'
    output:
        dea_fname = 'data/{source}/dea_result.csv',
        out_dir = directory('data/{source}/dea_misc/')
    script:
        'scripts/perform_dea.R'


rule analyze_dea:
    input:
        dea_fname = 'data/{source}/dea_result.csv'
    output:
        out_dir = directory('data/{source}/dea_stats/')
    log:
        notebook = 'notebooks/AnalyzeDEA.{source}.ipynb'
    notebook:
        'notebooks/AnalyzeDEA.ipynb'


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
            expression_mode \
            "$work_dir/{output.run_dir}" \
            "$work_dir/{output.result_fname}" \
            "$work_dir/{output.meta_fname}" \
            "$work_dir/{input.cell_fname}" \
            "$work_dir/{input.expr_fname}" \
            "$work_dir/{input.info_fname}" \
            "$work_dir/{input.dea_fname}"
        """


rule compare_tools:
    input:
        result_files = expand('runs/{tool}/{{source}}/enrichment_result.csv', tool=tool_list),
        meta_files = expand('runs/{tool}/{{source}}/meta.json', tool=tool_list)
    output:
        out_dir = directory('comparison/{source}/')
    log:
        notebook = 'notebooks/ToolComparison.{source}.ipynb'
    notebook:
        'notebooks/ToolComparison.ipynb'
