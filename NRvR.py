import os
import json
import numpy as np
import pandas as pd
import scipy.io
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
from scipy.stats import ks_2samp

# Set up files
MAT_FOLDER = "/Users/dineshchavan/Documents/BrainData/BRAIN_IMC_CellType"
SEG_FOLDER = "/Users/dineshchavan/Documents/BrainData/BRAIN_IMC_Segmentation"
OUTPUT_FOLDER = "/Users/dineshchavan/Downloads/metastasis_combined"
NRVR_JSON = os.path.join(OUTPUT_FOLDER, "nrvrV2_jsons")
COORDS_JSON_FOLDER = os.path.join(OUTPUT_FOLDER, "SlideCoordinates")

os.makedirs(NRVR_JSON, exist_ok=True)
os.makedirs(COORDS_JSON_FOLDER, exist_ok=True)

N_SHUFFLES = 5000
DBSCAN_EPS = 30
DBSCAN_MIN_SAMPLES = 3
MIN_CELL_COUNT = 20
ALPHA = 0.05

def ind2sub(array_shape, ind):
    rows = ind // array_shape[1]
    cols = ind % array_shape[1]
    return cols, rows

def load_mat_data(mat_path):
    slide_name = os.path.basename(mat_path).replace(".mat", "")
    seg_path = os.path.join(SEG_FOLDER, slide_name, "nuclei_multiscale.mat")

    try:
        seg_data = scipy.io.loadmat(seg_path)
        image_width = len(seg_data['nucleiImage'])
        image_size = (image_width, image_width)
    except:
        image_size = (1000, 1000)

    data = scipy.io.loadmat(mat_path)
    bounds = data['Boundaries'][0]
    types = data['cellTypes']
    coords, ctypes = [], []
    for i, arr in enumerate(bounds):
        lin = arr.flatten()
        xs, ys = ind2sub(image_size, lin)
        coords.append([np.mean(xs), np.mean(ys)])
        tarr = types[i]
        ctypes.append(tarr[0][0] if tarr.size and tarr[0].size else "Unknown")
    return np.array(coords), ctypes

def save_coords_json(image_id, coords, ctypes):
    cell_data = []
    for i, coord in enumerate(coords):
        cell_data.append({
            "phenotype": ctypes[i],
            "x": float(coord[0]),
            "y": float(coord[1])
        })
    out = {"slide_id": image_id, "cells": cell_data}
    with open(os.path.join(COORDS_JSON_FOLDER, f"{image_id}.json"), "w") as f:
        json.dump(out, f, indent=4)

# fencing metric
def dbscan_fraction_near_cancer(target_coords, cancer_coords, eps, min_samples):
    total = len(target_coords)
    if total < min_samples:
        return 0.0, 0, 0
    labels = DBSCAN(eps=eps, min_samples=min_samples).fit(target_coords).labels_
    unique = set(labels)
    total_clusters = sum(1 for l in unique if l != -1)
    valid_cells = valid_clusters = 0
    for lbl in unique:
        if lbl == -1:
            continue
        pts = target_coords[labels == lbl]
        if cancer_coords.size == 0:
            continue
        dists = np.linalg.norm(pts[:, None] - cancer_coords[None, :], axis=2)
        if np.any(dists <= eps):
            valid_cells += len(pts)
            valid_clusters += 1
    frac = valid_cells / total if total else 0.0
    return frac, valid_clusters, total_clusters


def run_nrv_analysis():
    mat_files = sorted([f for f in os.listdir(MAT_FOLDER) if f.endswith(".mat")])
    all_cell_types = set()
    for mat_file in mat_files:
        try:
            coords, ctypes = load_mat_data(os.path.join(MAT_FOLDER, mat_file))
            all_cell_types.update(ctypes)
        except:
            continue
    all_cell_types.discard("Cancer")
    all_cell_types = sorted(all_cell_types)

    print(f"\n Starting NRvR analysis on {len(all_cell_types)} cell types...\n")

    for cell_type in tqdm(all_cell_types, desc="Cell types", unit="type"):
        metastatic_values, nonmetastatic_values = [], []

        for mat_file in tqdm(mat_files, desc=f"Slides ({cell_type})", unit="file", leave=False):
            image_id = mat_file.replace(".mat", "")
            mat_path = os.path.join(MAT_FOLDER, mat_file)

            try:
                coords, ctypes = load_mat_data(mat_path)
                save_coords_json(image_id, coords, ctypes)
            except:
                continue

          ## Distinguish responder types through cell ids
            if image_id.startswith("BrM"):
                tumor_type = "Metastatic"
            elif image_id.startswith("Glioma"):
                tumor_type = "Non-Metastatic"
            else:
                tumor_type = "Unknown"

            tgt_mask = [ct == cell_type for ct in ctypes]
            can_mask = [ct == "Cancer" for ct in ctypes]
            other_mask = [(not t) and (not c) for t, c in zip(tgt_mask, can_mask)]

            tgt = coords[tgt_mask]
            can = coords[can_mask]
            other = coords[other_mask]

            if len(tgt) < MIN_CELL_COUNT or len(can) < MIN_CELL_COUNT or len(other) == 0:
                continue

            real_frac, _, _ = dbscan_fraction_near_cancer(tgt, can, DBSCAN_EPS, DBSCAN_MIN_SAMPLES)

            random_fracs = []
            for _ in range(N_SHUFFLES):
                idx = np.random.choice(len(other), size=len(tgt), replace=(len(tgt) > len(other)))
                sample = other[idx]
                shuffle_frac, _, _ = dbscan_fraction_near_cancer(sample, can, DBSCAN_EPS, DBSCAN_MIN_SAMPLES)
                random_fracs.append(shuffle_frac)

            avg_rand = np.mean(random_fracs)
            fencing_metric = (real_frac - avg_rand) / real_frac if real_frac > 0 else 0.0

            if tumor_type == "Metastatic":
                metastatic_values.append(fencing_metric)
            elif tumor_type == "Non-Metastatic":
                nonmetastatic_values.append(fencing_metric)

        if not metastatic_values and not nonmetastatic_values:
            continue

        if metastatic_values and nonmetastatic_values:
            ks_stat, ks_pval = ks_2samp(metastatic_values, nonmetastatic_values)
        else:
            ks_stat, ks_pval = None, None

        result = {
            "phenotype": cell_type,
            "n_valid_slides": len(metastatic_values) + len(nonmetastatic_values),
            "n_responders": len(metastatic_values),
            "n_nonresponders": len(nonmetastatic_values),
            "ks_statistic": ks_stat,
            "ks_pvalue": ks_pval,
            "responder_values": metastatic_values,
            "nonresponder_values": nonmetastatic_values
        }

        nrvr_path = os.path.join(NRVR_JSON, f"{cell_type.replace(' ', '_')}.json")
        with open(nrvr_path, "w") as f:
            json.dump(result, f, indent=4)

    print("\n Finished NRvR analysis!")

## main method
if __name__ == "__main__":
    run_nrv_analysis()
