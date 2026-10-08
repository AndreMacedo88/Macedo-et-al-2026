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

import os
import pathlib
import re
import pandas as pd


# anno_csv = pd.read_csv(config["anno_csv"], sep=",", index_col=0)
# ALL_IDS = anno_csv.iloc[:, 0]


configfile: os.path.join(workflow.basedir, "01-config_utrlevel.yaml")
config.update(_expand_env(dict(config)))


# def get_sample_ids(csv_file):
#     def check_matches(matches):
#         if len(matches) == 1:
#             batch_dirs.append(matches[0])
#         elif len(matches) == 0:
#             raise RuntimeError(f"No directory matched ID: {batch_id}")
#         else:
#             raise RuntimeError(f"More than one directory matched ID: {batch_id}")

#     def get_id_column_name(df):
#         id_column = df.columns[df.columns.str.lower().str.contains("id")]
#         if len(id_column) > 1:
#             raise RuntimeError(
#                 (
#                     f"Please make sure you have unique ID columns.\n",
#                     f"Currently the following columns are seen as ID columns: {id_column}",
#                 )
#             )
#         return id_column

#     df = pd.read_csv(csv_file, skip_blank_lines=True)
#     df = df.dropna(axis=0, how="all")
#     id_column = get_id_column_name(df)[0]

#     batch_dirs = []
#     for batch_id in df.Batch:
#         try:
#             matches = [
#                 dir_
#                 for dir_ in os.listdir(path=config["fastq_dir"])
#                 if batch_id in dir_
#             ]

#             check_matches(matches)
#         except RuntimeError as err:
#             batch_id = batch_id.replace(":", "_")
#             matches = [
#                 dir_
#                 for dir_ in os.listdir(path=config["fastq_dir"])
#                 if batch_id in dir_
#             ]

#             try:
#                 check_matches(matches)
#             except RuntimeError:
#                 raise err
#     return dict(zip(df[id_column], batch_dirs))


# def get_sample_paths(sample_ids):
#     sample_paths = dict()
#     for sample_id, batch_dir in sample_ids.items():
#         try:
#             sample_paths[sample_id] = [
#                 str(
#                     next(
#                         pathlib.Path(
#                             os.path.join(config["fastq_dir"], batch_dir)
#                         ).rglob(f"{sample_id}_R1.fastq.gz")
#                     )
#                 )
#             ]
#             sample_paths[sample_id].append(
#                 str(
#                     next(
#                         pathlib.Path(
#                             os.path.join(config["fastq_dir"], batch_dir)
#                         ).rglob(f"{sample_id}_R2.fastq.gz")
#                     )
#                 )
#             )
#         except StopIteration:
#             print(f"Read 1 file for sample {sample_id} of batch {batch_dir} not found.")
#     return sample_paths


# def get_input_dirs_r1(wildcards):
#     return IDS_AND_PATHS.get(wildcards.sample_id)[0]


# def get_input_dirs_r2(wildcards):
#     return IDS_AND_PATHS.get(wildcards.sample_id)[1]


# sample_ids = get_sample_ids(config["anno_csv"])
# IDS_AND_PATHS = get_sample_paths(sample_ids)
# ALL_IDS = IDS_AND_PATHS.keys()


def get_data(folders_path):
    folder_dict = dict()
    file_list = []
    for path in pathlib.Path(folders_path).rglob("*R1.fastq.gz"):
        id_name = re.match(f"{folders_path}(.*\d)_R\d", str(path))
        try:
            folder_name = id_name.group(1)
            file_name = folder_name.split("/")[2]
            file_list.append(file_name)
            folder_dict[file_name] = folder_name
        except AttributeError:
            continue
    return file_list, folder_dict


ALL_IDS, ALL_FOLDERS = get_data(f"{config['fastq_dir']}")


def get_input_dirs_r1(wildcards):
    input_dir = "".join(
        [config["fastq_dir"], ALL_FOLDERS.get(wildcards.sample_id), "_R1.fastq.gz"]
    )
    return input_dir


def get_input_dirs_r2(wildcards):
    input_dir = "".join(
        [config["fastq_dir"], ALL_FOLDERS.get(wildcards.sample_id), "_R2.fastq.gz"]
    )
    return input_dir


rule all:
    input:
        f"{config['results_dir']}star/mapping.csv",
        f"{config['results_dir']}{config['results_base_name']}_bw12_above10cts_apa_named.h5ad",
        f"{config['results_dir']}{config['results_base_name']}_means.h5ad",


rule star_map:
    input:
        r1=get_input_dirs_r1,
        r2=get_input_dirs_r2,
        index=f"{config['index_star']}",
    output:
        "{dir_results}star/{sample_id}/Aligned.sortedByCoord.out.bam",
        "{dir_results}star/{sample_id}/Log.final.out",
        "{dir_results}star/{sample_id}/Log.out",
        "{dir_results}star/{sample_id}/Log.progress.out",
        "{dir_results}star/{sample_id}/SJ.out.tab",
    threads: 26
    params:
        basename="{dir_results}star/{sample_id}/",
    conda:
        "01-mapping_star.yaml"
    # As originally run, the strings for --outFilterMismatchNmax, --outFilterMismatchNoverReadLmax and
    # --chimOutType lacked separating spaces, so STAR applied only --outFilterMismatchNmax 15 and kept the
    # defaults for the other two (checked in the STAR Log.out files). The command below states the
    # parameters that were in effect.
    shell:
        "STAR "
        "--runThreadN {threads} "
        "--genomeDir {input.index} "
        "--outFileNamePrefix {params.basename} "
        "--outSAMtype BAM SortedByCoordinate "
        "--outSAMattributes All "
        "--outSAMunmapped Within "
        "--bamRemoveDuplicatesType UniqueIdentical "
        "--outFilterMatchNminOverLread 0.2 "
        "--outFilterScoreMinOverLread 0.3 "
        "--outFilterMismatchNmax 15 "
        "--readFilesCommand zcat "
        "--readFilesIn {input.r1} {input.r2} "


rule samtools_index:
    input:
        "{dir_results}star/{sample_id}/Aligned.sortedByCoord.out.bam",
    output:
        "{dir_results}star/{sample_id}/Aligned.sortedByCoord.out.bam.bai",
    conda:
        "01-mapping_star.yaml"
    shell:
        "samtools index {input}"


rule symlink_bams:
    input:
        f"{config['results_dir']}star/{{sample_id}}/Aligned.sortedByCoord.out.bam",
    output:
        f"{config['results_dir']}star/bam_lns/{{sample_id}}.sorted.bam",
    shell:
        "ln -s {input} {output}"


rule symlink_bais:
    input:
        f"{config['results_dir']}star/{{sample_id}}/Aligned.sortedByCoord.out.bam.bai",
    output:
        f"{config['results_dir']}star/bam_lns/{{sample_id}}.sorted.bam.bai",
    shell:
        "ln -s {input} {output}"


rule mapping_rate:
    input:
        star=expand(
            "{dir_results}star/{sample_id}/Log.final.out",
            dir_results=config["results_dir"],
            sample_id=ALL_IDS,
        ),
    params:
        ids=ALL_IDS,
    output:
        csv=f"{config['results_dir']}star/mapping.csv",
    run:
        import pandas as pd

        data = pd.DataFrame({"id": params["ids"], "star_mapfile": input.star})
        data["data"] = list(
            map(partial(pd.read_csv, sep="\t", header=None), data["star_mapfile"])
        )
        data["data"] = list(map(lambda x: x[~x[1].isna()], data.data))
        data["data"] = list(
            map(
                lambda x: {k.strip(" |"): v.strip("%") for (k, v) in zip(x[0], x[1])},
                data.data,
            )
        )
        data = data.join(pd.DataFrame(list(data.data)))
        data.drop(columns="data").to_csv(output.csv, index=False)


rule utr_mapping_raw:
    input:
        bams=expand(
            f"{config['results_dir']}star/bam_lns/{{sample_id}}.sorted.bam",
            sample_id=ALL_IDS,
        ),
        bais=expand(
            f"{config['results_dir']}star/bam_lns/{{sample_id}}.sorted.bam.bai",
            sample_id=ALL_IDS,
        ),
    output:
        "{dir_results}/{results_base_name}_raw.h5ad",
    params:
        gtf=config["human_gtf"],
        anno=config["anno_csv"],
    threads: 26
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling run {params.gtf} -a {input.bams} -f three_prime_utr -s utr -p smartseq -v 3 -m union -r -A {params.anno} -L True -i -o {output}
        """


# sample_range = range(1, 20)  # from 1 until (sample number + 1)


# rule utr_mapping_raw_callpeaks:
#     input:
#         counts=expand(
#             f"intermediate/UTR_sizes_counts_{{count}}.parquet", count=sample_range
#         ),
#     output:
#         "{dir_results}/{results_base_name}_raw.h5ad",
#     params:
#         utr_regions="intermediate/utr_reference_regions.parquet",
#         anno=config["anno_csv"],
#     threads: 26
#     conda:
#         "apa-utrcalling"
#     shell:
#         """
#         utrcalling call-peaks {params.utr_regions} -c {input.counts} -s utr -p smartseq -v 3 -m union -r -A {params.anno} -L True -i -o {output}
#         """


# Use the next rule to skip the first part of the UTR calling if some intermediate steps failed
# rule utr_mapping_raw_callpeaks_merged:
#     input:
#         counts="intermediate/UTR_sizes_counts_merged_union.parquet",
#     output:
#         "{dir_results}/{results_base_name}_raw.h5ad",
#     params:
#         utr_regions="intermediate/utr_reference_regions.parquet",
#         anno=config["anno_csv"],
#     threads: 26
#     conda:
#         "apa-utrcalling"
#     shell:
#         """
#         utrcalling call-peaks {params.utr_regions} -c {input.counts} -s utr -p smartseq -v 3 -m union -r -A {params.anno} -L True -i -o {output}
#         """


rule calculate_mean_utr_sizes:
    input:
        "{dir_results}{results_base_name}_raw.h5ad",
    output:
        "{dir_results}{results_base_name}_means.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling compute-mean-utr-sizes {input} {output} -u utr_name
        """


rule call_peaks:
    input:
        "{dir_results}{results_base_name}_raw.h5ad",
    output:
        "{dir_results}{results_base_name}_bw12.h5ad",
    params:
        out="{dir_results}{results_base_name}",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling call-peaks-different-bandwidths {input} {params.out} -b 12
        """


rule trim_peaks:
    input:
        "{dir_results}{results_base_name}_bw12.h5ad",
    output:
        "{dir_results}{results_base_name}_bw12_above10cts.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling remove-peaks-below-threshold {input} {output} -t 0.05 -m 10
        """


rule remove_non_apa_utrs:
    input:
        "{dir_results}{results_base_name}_bw12_above10cts.h5ad",
    output:
        "{dir_results}{results_base_name}_bw12_above10cts_apa.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling select-alternative-transcription-events {input} {output}
        """


rule rename_peaks_to_gene_name:
    input:
        "{dir_results}{results_base_name}_bw12_above10cts_apa.h5ad",
    output:
        "{dir_results}{results_base_name}_bw12_above10cts_apa_named.h5ad",
    conda:
        "apa-utrcalling"
    shell:
        """
        utrcalling generate-names {input} {output} -u utr_name -g gene_name
        """
