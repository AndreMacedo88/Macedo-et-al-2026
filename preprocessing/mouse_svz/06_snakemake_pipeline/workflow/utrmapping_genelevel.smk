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

configfile: os.path.join(workflow.basedir, "..", "config", "config_genelevel.yaml")
config.update(_expand_env(dict(config)))

ALL_VERSIONS = ["2", "2", "3", "3", "3", "3"]
SUFFIX = ["1", "2", "3", "4", "5", "6"]
ALL_SAMPLES = [
    "AS-196777-LR-30151",
    "AS-196779-LR-30152",
    "wildtype1",
    "wildtype2",
    "wildtype3",
    "wildtype4",
]


def get_bams(wildcards):
    return expand(
        f"{config['results_dir']}data/{{sample_id}}.sorted.bam",
        sample_id=ALL_SAMPLES[SUFFIX.index(wildcards.suffix)],
    )


def get_bais(wildcards):
    return expand(
        f"{config['results_dir']}data/{{sample_id}}.sorted.bam.bai",
        sample_id=ALL_SAMPLES[SUFFIX.index(wildcards.suffix)],
    )


def get_anno(wildcards):
    return config[f"anno_csv_v{wildcards.version}_{wildcards.suffix}"]


rule all:
    input:
        f"{config['results_dir']}utrcalling/{config['results_base_name']}_bw12_above10cts_apa_named_percentages.h5ad",


rule utr_mapping_raw:
    input:
        bams=get_bams,
        bais=get_bais,
    output:
        out = "{dir_results}utrcalling/{results_base_name}_{suffix}_v{version}.h5ad",
        peaks = "{dir_results}utrcalling/{results_base_name}_{suffix}_v{version}_kde.h5ad",
    params:
        gtf=config["mouse_gtf"],
        anno=get_anno,
    threads: 26
    conda:
        "utrcalling"
    shell:
        """
        utrcalling run {params.gtf} -a {input.bams} -f three_prime_utr -s gene -p 10x -v {wildcards.version} -m union -r -A {params.anno} -L False -i -o {output.out}
        """


def merge_utr_mappings_input(wildcards):
    input_merge = [
        f"{wildcards.dir_results}utrcalling/{wildcards.results_base_name}_{suffix}_v{version}.h5ad"
        for suffix, version in zip(SUFFIX, ALL_VERSIONS)
    ]
    return input_merge


rule merge_utr_mappings:
    input:
        merge_utr_mappings_input,
    output:
        "{dir_results}utrcalling/{results_base_name}_raw.h5ad",
    conda:
        "downstr"
    shell:
        "python {workflow.basedir}/merge_utr_mappings.py {input} -o {output}"


rule call_peaks:
    input:
        "{dir_results}utrcalling/{results_base_name}_raw.h5ad",
    output:
        "{dir_results}utrcalling/{results_base_name}_bw12.h5ad",
    params:
        out="{dir_results}utrcalling/{results_base_name}",
    conda:
        "utrcalling"
    shell:
        """
        utrcalling call-peaks-different-bandwidths {input} {params.out} -b 12
        """

rule trim_peaks:
    input:
        "{dir_results}utrcalling/{results_base_name}_bw12.h5ad",
    output:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts.h5ad",
    conda:
        "utrcalling"
    shell:
        """
        utrcalling remove-peaks-below-threshold {input} {output} -t 0.05 -m 10
        """


rule remove_non_apa_utrs:
    input:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts.h5ad",
    output:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts_apa.h5ad",
    conda:
        "utrcalling"
    shell:
        """
        utrcalling select-alternative-transcription-events {input} {output}
        """


rule rename_peaks_to_gene_name:
    input:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts_apa.h5ad",
    output:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts_apa_named.h5ad",
    conda:
        "utrcalling"
    shell:
        """
        utrcalling generate-names {input} {output} -u utr_name -g gene_name
        """

rule counts_to_percentages:
    input:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts_apa_named.h5ad",
    output:
        "{dir_results}utrcalling/{results_base_name}_bw12_above10cts_apa_named_percentages.h5ad",
    conda:
        "utrcalling"
    shell:
        """
        utrcalling counts-to-percentages {input} {output} -u utr_name
        """
