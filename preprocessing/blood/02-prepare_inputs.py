import os
import re
import scanpy as sc

DATA_ROOT = os.environ["M26_DATA_ROOT"]  # data directory (see README)

BASEPATH = f"{DATA_ROOT}/3utr/Blood/"
PREPROCESSINGPATH = {
    "all_till_06Jul": os.path.join(BASEPATH, "all_till_06Jul/data"),
    "20230725": os.path.join(BASEPATH, "20230725"),
    "20230816": os.path.join(BASEPATH, "20230816"),
    "20240724": os.path.join(BASEPATH, "20240724"),
}
RESULTSPATH = os.path.join(BASEPATH, "all_till_20240724/utrcalling/")
PROCESSINGPATH = os.path.join(BASEPATH, "all_till_20240724/processing/")

bam_paths = []
for dir_name, preprocess_path in PREPROCESSINGPATH.items():
    pattern = re.compile(os.path.join(preprocess_path, "star/(.*)/Aligned.sortedByCoord.out.bam"))

    for root, dirs, files in os.walk(os.path.join(preprocess_path, "star"), followlinks=True):
        for file in files:
            file_path = os.path.join(root, file)
            match = pattern.fullmatch(file_path)
            if match:
                bam_paths.append(match.string)


annotations_file_full = sc.read_h5ad(
    os.path.join(
        RESULTSPATH,
        "UTR_count_table_UTRlevel_peaks_handmerged_bw12_100cts_bw12_trimmed_APA_proportions_fullobs.h5ad",
    )
).obs

annotations_file = annotations_file_full.drop("merge_batch", axis=1)
annotations_file.Batch = annotations_file.Batch.astype(str)
annotations_file.read_barcodes = annotations_file.read_barcodes.astype(str)
annotations_file["read_bcs_old"] = annotations_file.read_barcodes.copy()
annotations_file.read_barcodes = annotations_file.read_barcodes.str.split("-").str[0]

sample_ids = set(annotations_file["sample_id"])
bam_paths = [
    bam_path for bam_path in bam_paths if any(sample_id in bam_path for sample_id in sample_ids)
]

with open(os.path.join(PROCESSINGPATH, "bam_paths.txt"), "w") as f:
    f.write("\n".join(bam_paths))

annotations_file.to_csv(os.path.join(PROCESSINGPATH, "annotations.csv"))
