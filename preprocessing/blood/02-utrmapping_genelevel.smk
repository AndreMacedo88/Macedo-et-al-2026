import os

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)


def _expand_env(value):
    """Expand ${M26_DATA_ROOT} in config values."""
    if isinstance(value, str):
        return os.path.expandvars(value)
    if isinstance(value, dict):
        return {k: _expand_env(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_expand_env(v) for v in value]
    return value

configfile: os.path.join(workflow.basedir, "02-config_genelevel.yaml")
config.update(_expand_env(dict(config)))


rule all:
    input:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa_named_percentages.h5ad",


rule prepare_inputs:
    output:
        bam_list=f"{config['process_dir']}/bam_paths.txt",
        anno_file=f"{config['process_dir']}/annotations.csv",
    conda:
        "apa-python"
    script:
        "02-prepare_inputs.py"


rule utr_mapping_raw:
    input:
        bam_list=rules.prepare_inputs.output.bam_list,
        anno_file=rules.prepare_inputs.output.anno_file,
    output:
        "{dir_results}/{results_base_name}_raw.h5ad",
    params:
        gtf=config["human_gtf"],
    threads: 26
    conda:
        "apa-utrcalling"
    shell:
        """
        bams=$(cat {input.bam_list} | tr '\n' ' ')
        utrcalling run {params.gtf} -a $bams -f three_prime_utr -s gene -p smartseq -v 3 -m union -r -A {input.anno_file} -L True -i -o {output}
        """


rule call_peaks:
    input:
        config["raw_utr_lengths"],
    output:
        f"{config['results_dir']}{config['results_base_name']}_bw12.h5ad",
    params:
        out=f"{config['results_dir']}{config['results_base_name']}",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling call-peaks-different-bandwidths {input} {params.out} -b 12
        """


rule trim_peaks:
    input:
        f"{config['results_dir']}{config['results_base_name']}_bw12.h5ad",
    output:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling remove-peaks-below-threshold {input} {output} -t 0.05 -m 10
        """


rule remove_non_apa_utrs:
    input:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts.h5ad",
    output:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling select-alternative-transcription-events {input} {output}
        """


rule rename_peaks_to_gene_name:
    input:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa.h5ad",
    output:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa_named.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling generate-names {input} {output} -u utr_name -g gene_name
        """

rule counts_to_percentages:
    input:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa_named.h5ad",
    output:
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa_named_percentages.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling counts-to-percentages {input} {output} -u utr_name
        """