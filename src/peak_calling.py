"""
Functions for clustering data points into peaks and computing PAS distances.
"""
# import general tools
import os
from urllib.error import HTTPError
import time
import regex
# import numerical methods
import numpy as np
import math
from sklearn.neighbors import KernelDensity
from scipy.signal import argrelextrema
# import data processing tools
import pandas as pd
import anndata as ad
# import biological tools
from Bio import Entrez, SeqIO
import pybedtools
# import quality of life tools
import joblib
from tqdm import tqdm


def peak_call_anndata_parallel(adata, method, **kwargs):

    def peak_call(adata, method, utr):
        adata_utr = adata[:, adata.var.utr_name == utr]
        gene = adata_utr.var.gene_name.iloc[0]
        utr_name_parts = adata_utr.var.utr_name.iloc[0].split("_")
        if len(utr_name_parts) > 1:
            utr_number = f"UTR{utr_name_parts[1]}_"
        else:
            utr_number = ""
        utr_var = adata_utr.var.iloc[0]

        # generate a vector of one UTR size per UMI count
        umi_utr_sizes = _utr_sizes_vector_from_read_counts(
            adata_utr.X, adata_utr.var.utr_size)

        if str.lower(method) == "kde":
            peaks = _call_peaks_using_kernel_density_estimator(
                umi_utr_sizes.values, bandwidth=kwargs.get("bandwidth", 1)
            )
            peak_sizes = pd.Series(data=peaks, index=umi_utr_sizes.index)
            peak_sizes = peak_sizes[~peak_sizes.index.duplicated()]
            peak_sizes_unique = peak_sizes.unique()

        # create new matrix with the sum of UMIs for each new peak
        mtx_new_t = adata_utr.to_df().T
        mtx_new_t["sizes"] = peak_sizes.tolist()
        mtx_new = mtx_new_t.groupby("sizes").sum().T

        # create new polyA names
        umi_utr_sizes_new = mtx_new.columns
        polya_names = [f"{gene}_{utr_number}pA{x}" for x in range(
            1, len(umi_utr_sizes_new) + 1)]

        # generate new utr_var
        utr_var_new = pd.DataFrame([utr_var] * len(umi_utr_sizes_new))
        utr_var_new.index = polya_names
        utr_var_new["utr_size"] = peak_sizes_unique

        return mtx_new.to_numpy(), utr_var_new

    if kwargs.get("bandwidth", 1) == 0:
        print("Using Bandwidth = 0 means no changes were made")
        return adata

    utrs = adata.var.utr_name.unique()

    parallel_results = joblib.Parallel(n_jobs=-2)(
        joblib.delayed(peak_call)(adata=adata, method=method, utr=utr) for utr in tqdm(utrs)
    )

    mtx_list, var_list = list(zip(*parallel_results))

    var_df = pd.concat(var_list, axis=0).drop("pA", axis=1, errors="ignore")
    new_adata = ad.AnnData(X=np.hstack(mtx_list), obs=adata.obs, var=var_df)
    new_adata = new_adata[new_adata.obs.sort_index().index, :]
    return new_adata


def peak_call_anndata(adata, method, **kwargs):
    if kwargs.get("bandwidth", 1) == 0:
        print("Using Bandwidth = 0 means no changes were made")
        return adata

    utrs = adata.var.utr_name.unique()

    mtx_list = []
    var_list = []

    # select the UTR
    for utr in tqdm(utrs, desc="Calling peaks for UTR: "):
        adata_utr = adata[:, adata.var.utr_name == utr]
        gene = adata_utr.var.gene_name.iloc[0]
        utr_name_parts = adata_utr.var.utr_name.iloc[0].split("_")
        if len(utr_name_parts) > 1:
            utr_number = f"UTR{utr_name_parts[1]}_"
        else:
            utr_number = ""
        utr_var = adata_utr.var.iloc[0]

        # generate a vector of one UTR size per UMI count
        umi_utr_sizes = _utr_sizes_vector_from_read_counts(
            adata_utr.X, adata_utr.var.utr_size)

        if str.lower(method) == "kde":
            peaks = _call_peaks_using_kernel_density_estimator(
                umi_utr_sizes.values, bandwidth=kwargs.get("bandwidth", 1)
            )
            peak_sizes = pd.Series(data=peaks, index=umi_utr_sizes.index)
            peak_sizes = peak_sizes[~peak_sizes.index.duplicated()]
            peak_sizes_unique = peak_sizes.unique()

        # create new matrix with the sum of UMIs for each new peak
        mtx_new_t = adata_utr.to_df().T
        mtx_new_t["sizes"] = peak_sizes.tolist()
        mtx_new = mtx_new_t.groupby("sizes").sum().T
        mtx_list.append(mtx_new.to_numpy())

        # create new polyA names
        umi_utr_sizes_new = mtx_new.columns
        polya_names = [f"{gene}_{utr_number}pA{x}" for x in range(
            1, len(umi_utr_sizes_new) + 1)]

        # generate new utr_var
        utr_var_new = pd.DataFrame([utr_var] * len(umi_utr_sizes_new))
        utr_var_new.index = polya_names
        utr_var_new["utr_size"] = peak_sizes_unique
        var_list.append(utr_var_new)

    var_df = pd.concat(var_list, axis=0).drop("pA", axis=1, errors="ignore")
    new_adata = ad.AnnData(X=np.hstack(mtx_list), obs=adata.obs, var=var_df)
    new_adata = new_adata[new_adata.obs.sort_index().index, :]
    return new_adata


def remove_peaks_below_count(adata, threshold=10):
    """
    Removes the UTR location peaks that have less than `threshold` counts in
    the sum of all samples. This effectively curates the peaks, selecting the most common
    and/or more prone to be true peaks.

    Args:
        adata (ad.AnnData): The AnnData object to filter.
        threshold (int, optional): The minimum counts number each peak must have.
            Defaults to 10.

    Returns:
        ad.AnnData: A filtered AnnData object with UTR peaks counts above the threshold.
    """
    col_sums = np.asarray(adata.X.sum(axis=0)).ravel()
    mask = col_sums >= threshold
    return adata[:, mask].copy()


def estimate_best_bandwidth(adata, method="kde", genome=None, hamming=1,
                            bandwidths=(0.5, 1, 1.5, 2, 2.5, 3),
                            return_distances_to_hex=False):

    if isinstance(bandwidths, int):
        bandwidths = [bandwidths]

    peak_data = joblib.Parallel(n_jobs=-2)(
        joblib.delayed(peak_call_anndata)(adata=adata, method=method, bandwidth=bandwidth) for bandwidth in bandwidths
    )

    estimation_results = joblib.Parallel(n_jobs=-2)(
        joblib.delayed(_score_peak_assignment)(adata=adata_bandwidth, genome=genome,
                                               hamming=hamming,
                                               return_distances_to_hex=return_distances_to_hex) for adata_bandwidth in peak_data
    )
    return estimation_results


def pas_distance(adata, genome, hamming=1):
    """
    Calculates the distance of each read to the nearest polyadenilation hexamer site
    AATAAA.

    Args:
        adata: an anndata file containing UMI counts per UTR, with UTR sizes
            in .var.utr_size and UTR names in .var.utr_name.
        genome: a fasta file of a genome assembly (pybedtools.BedTool).
        hamming (int, optional): value for the maximum degeneration of the PAS.
            Defaults to 1.
    """

    def get_distances_per_utr(utr):
        adata_utr = adata[:, adata.var.utr_name == utr]

        pas_end, _, _ = _find_pas(adata_utr, genome=genome, hamming=hamming)

        if not pas_end:
            return None

        # generate a vector of one UTR size per UMI count
        umi_utr_sizes = _utr_sizes_vector_from_read_counts(
            adata_utr.X, adata_utr.var.utr_size)

        distances = _calculate_umi_distance_to_pas(
            umi_utr_sizes=np.asarray(umi_utr_sizes), pas_end=pas_end)
        return np.asarray(distances).ravel()

    utrs = adata.var.utr_name.unique()

    distances_list = joblib.Parallel(n_jobs=-2)(
        joblib.delayed(get_distances_per_utr)(utr=utr) for utr in tqdm(utrs)
    )

    # remove None from the list
    distances_list = [i for i in distances_list if i is not None]

    distances_array = np.concatenate(distances_list)

    return distances_array


def _score_peak_assignment(adata, genome=None, hamming=1, return_distances_to_hex=False):

    def get_distances(utr):
        adata_utr = adata[:, adata.var.utr_name == utr]
        pas_end, _, _ = _find_pas(adata_utr, genome=genome, hamming=hamming)

        if not pas_end:
            return None

        # generate a vector of one UTR size per UMI count
        umi_utr_sizes = _utr_sizes_vector_from_read_counts(
            adata_utr.X, adata_utr.var.utr_size)

        distances_arr = _calculate_umi_distance_to_pas(
            umi_utr_sizes=np.asarray(umi_utr_sizes), pas_end=pas_end)
        distances_arr = np.asarray(distances_arr).ravel()

        total_umis = len(distances_arr)
        threshold = int(math.ceil(total_umis * 0.01))

        # Group by feature: count UMIs per UTR feature, keep features above threshold
        feature_labels = np.asarray(umi_utr_sizes.index) if isinstance(umi_utr_sizes, pd.Series) else np.zeros(total_umis)
        _, label_indices, label_counts = np.unique(
            feature_labels, return_inverse=True, return_counts=True)
        umi_feature_count = label_counts[label_indices]
        valid_mask = umi_feature_count > threshold
        distances_filtered = distances_arr[valid_mask]

        in_range = np.sum((distances_filtered >= 10) & (distances_filtered <= 30))

        return [in_range / total_umis if total_umis > 0 else 0, distances_arr]

    utrs = adata.var.utr_name.unique()

    results_per_utr = joblib.Parallel(n_jobs=-2)(
        joblib.delayed(get_distances)(utr=utr) for utr in tqdm(utrs, desc="Scoring peaks for UTR: ")
    )

    scores = []
    distances_list = []

    for i in results_per_utr:
        if i is None:
            continue
        scores.append(i[0])
        distances_list.append(i[1])

    score = sum(scores) / len(utrs) if len(utrs) > 0 else 0
    if return_distances_to_hex:
        distances_array = np.concatenate(distances_list)
        return (score, distances_array)
    else:
        return score


def _find_pas(adata, genome, hamming=1):
    seq_chrom = adata.var.chrom.iloc[0]
    seq_start = int(adata.var.start.iloc[0]) - 1
    seq_stop = adata.var.end.iloc[0]
    seq_strand = adata.var.strand.iloc[0]

    query = pybedtools.BedTool(
        f"{seq_chrom} {seq_start} {seq_stop} three_prime_utr 0 {seq_strand}", from_string=True)
    query = query.sequence(fi=genome, s=True)

    with open(query.seqfn) as q:
        next(q)
        # The next line gives us the sequence from the genome, reverse-complemented if necessary because of negative strand
        utr_seq = q.read().strip()

    # Compile regex once
    if hamming > 0:
        pattern = regex.compile(fr"(AATAAA){{s<={hamming}}}")
    elif hamming == 0:
        pattern = regex.compile(r"AATAAA")
    else:
        raise ValueError(f"Hamming distance of {hamming} not supported")

    pas_seq = []
    pas_end = []
    for match in pattern.finditer(utr_seq, overlapped=True):
        pas_end.append(match.end())
        pas_seq.append(match.group(0))

    return (pas_end, pas_seq, seq_strand)


def _find_pas_online(adata, hamming=1):
    chr_dict = {str(i): f"NC_{i:06d}" for i in range(1, 23)}
    chr_dict["X"] = "NC_000023"
    chr_dict["Y"] = "NC_000024"

    seq_chrom = adata.var.chrom.iloc[0]
    seq_chrom = chr_dict.get(seq_chrom, seq_chrom)
    seq_start = adata.var.start.iloc[0]
    seq_stop = adata.var.end.iloc[0]
    seq_strand = adata.var.strand.iloc[0]
    if seq_strand == "+":
        seq_strand = 1
    else:
        seq_strand = 2

    Entrez.email = 'anonymous@gmail.com'

    httperror_counter = 0
    while True:
        try:
            handle = Entrez.efetch(db="nucleotide",
                                   id=seq_chrom,
                                   rettype="fasta",
                                   strand=seq_strand,
                                   seq_start=seq_start,
                                   seq_stop=seq_stop
                                   )
            break
        except HTTPError as e:
            # if multiprocessing this function, there can be too many requests
            httperror_counter += 1
            if httperror_counter > 100:
                raise e
            else:
                time.sleep(2)

    record = SeqIO.read(handle, "fasta")
    handle.close()
    utr_seq = str(record.seq)

    # Compile regex once
    if hamming > 0:
        pattern = regex.compile(fr"(AATAAA){{s<={hamming}}}")
    elif hamming == 0:
        pattern = regex.compile(r"AATAAA")
    else:
        raise ValueError(f"Hamming distance of {hamming} not supported")

    pas_seq = []
    pas_end = []
    for match in pattern.finditer(utr_seq, overlapped=True):
        pas_end.append(match.end())
        pas_seq.append(match.group(0))

    return (pas_end, pas_seq, seq_strand)


def _calculate_umi_distance_to_pas(umi_utr_sizes, pas_end):
    """Calculate distance from each UMI to the nearest PAS site.

    Since pybedtools sequence extraction is strand-aware (s=True), the
    sequence is already oriented 5'->3' regardless of genomic strand.
    Distances are always: utr_size - pas_position.

    Args:
        umi_utr_sizes: 1D array-like of UTR sizes per UMI.
        pas_end: list of PAS end positions within the UTR sequence.

    Returns:
        numpy array of distances (one per UMI).
    """
    sizes = np.asarray(umi_utr_sizes).ravel()

    if len(pas_end) == 1:
        return sizes - pas_end[0]

    # Vectorized multi-PAS: compute distance to each PAS, pick closest
    pas_arr = np.array(pas_end)
    # Shape: (n_umis, n_pas)
    distances_matrix = sizes[:, np.newaxis] - pas_arr[np.newaxis, :]
    closest_idx = np.abs(distances_matrix).argmin(axis=1)
    return distances_matrix[np.arange(len(sizes)), closest_idx]


def _utr_sizes_vector_from_read_counts(counts_arr, sizes_arr):
    """Create a vector of UTR sizes, each repeated by its total UMI count.

    Handles both dense and sparse count matrices. Returns a pandas Series
    (preserving index from sizes_arr) when sizes_arr is a Series.

    Args:
        counts_arr: UMI count matrix (n_cells x n_features), dense or sparse.
        sizes_arr: pandas Series or array of UTR sizes per feature.

    Returns:
        pandas Series (if sizes_arr is Series) or numpy array of repeated sizes.
    """
    # .sum(axis=0) on sparse returns np.matrix (2D); flatten to 1D array
    counts_sum = np.asarray(counts_arr.sum(axis=0)).ravel()
    counts_sum = np.rint(counts_sum).astype(int)
    if isinstance(sizes_arr, pd.Series):
        return sizes_arr.repeat(counts_sum)
    return np.repeat(np.asarray(sizes_arr), counts_sum)


def _call_peaks_using_kernel_density_estimator(utr_sizes_vector, bandwidth=1):
    """
    Fit a kernel density model to predict the peak's distribution density function
    and assign each UMI to the nearest peak.

    Args:
        utr_sizes_vector (array-like): vector with the UTR sizes to call peaks.
        bandwidth (float): KDE bandwidth parameter.

    Returns:
        np.ndarray: The vector with the peak assignments.
    """
    utr_sizes_vector = np.asarray(utr_sizes_vector).ravel()
    n = len(utr_sizes_vector)

    if n <= 3:
        return np.full(n, round(np.median(utr_sizes_vector)))

    kde = KernelDensity(kernel='gaussian', bandwidth=bandwidth).fit(
        utr_sizes_vector.reshape(-1, 1))

    # calculate an appropriate window to check the model
    minimum_size, maximum_size = utr_sizes_vector.min(), utr_sizes_vector.max()
    padding = (maximum_size - minimum_size) * 0.1
    window = np.arange(round(minimum_size - padding), round(maximum_size + padding) + 1)

    # get the model predictions for that window
    evaluation = kde.score_samples(window.reshape(-1, 1))
    minima_pos = window[argrelextrema(evaluation, np.less)[0]]
    maxima_pos = window[argrelextrema(evaluation, np.greater)[0]]

    if len(maxima_pos) == 0:
        return np.full(n, round(np.median(utr_sizes_vector)))

    # Assign each UMI to a bin defined by the KDE minima, then map to peak center
    bin_indices = np.digitize(utr_sizes_vector, minima_pos)
    bin_indices = np.clip(bin_indices, 0, len(maxima_pos) - 1)
    peaks = np.round(maxima_pos[bin_indices]).astype(float)

    return peaks
