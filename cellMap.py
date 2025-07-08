import os
import scipy.io
import numpy as np
import matplotlib.pyplot as plt


MAT_FOLDER    = "BrainData/BRAIN_IMC_CellType"
OUTPUT_FOLDER = "CellTypePlots"
PLOT_FOLDER   = os.path.join(OUTPUT_FOLDER, "GliomaAllSlidesPlot")

os.makedirs(PLOT_FOLDER, exist_ok=True)

IMAGE_SIZE = (1000, 1000)


def ind2sub(array_shape, ind):
    rows = ind // array_shape[1]
    cols = ind % array_shape[1]
    return cols, rows

def load_mat_data(mat_path):
    d = scipy.io.loadmat(mat_path)
    bounds = d['Boundaries'][0]
    types  = d['cellTypes']
    coords = []
    ctypes = []
    for i, arr in enumerate(bounds):
        lin = arr.flatten()
        x, y = ind2sub(IMAGE_SIZE, lin)
        coords.append([np.mean(x), np.mean(y)])
        tarr = types[i]
        if tarr.size and tarr[0].size and len(tarr[0][0]):
            ctypes.append(tarr[0][0])
        else:
            ctypes.append("Unknown")
    return np.array(coords), ctypes

# 1) Build a global list of all cell‑types across every slide
all_types = set()
for fn in os.listdir(MAT_FOLDER):
    if not fn.endswith(".mat"):
        continue
    coords, ctypes = load_mat_data(os.path.join(MAT_FOLDER, fn))
    all_types.update(ctypes)


palette = {
    "Tc":               "blue",
    "Cancer":           "red",
    "Endothelial cell": "purple",
    "Unknown":          "gray"
}

others = sorted(t for t in all_types if t not in palette)
cmap   = plt.cm.get_cmap("tab20", len(others))
for idx, ct in enumerate(others):
    palette[ct] = cmap(idx)

# 3) Loop through each slide, plot & save
for fn in sorted(os.listdir(MAT_FOLDER)):
    if not fn.endswith(".mat"):
        continue

    slide_id = os.path.splitext(fn)[0]
    mat_path  = os.path.join(MAT_FOLDER, fn)

    coords, ctypes = load_mat_data(mat_path)

    # Group indices by cell type
    groups = {}
    for i, ct in enumerate(ctypes):
        groups.setdefault(ct, []).append(i)

    # Make the plot
    fig, ax = plt.subplots(figsize=(8,8))
    for ct, idxs in groups.items():
        pts = coords[idxs]
        ax.scatter(
            pts[:,0], pts[:,1],
            c=palette.get(ct, "black"),
            s=15,
            label=ct,
            alpha=0.8,
            marker='x' if ct=="Cancer" else 'o'
        )

    ax.set_title(f"All cell types — {slide_id}")
    ax.set_xlabel("X (μm)")
    ax.set_ylabel("Y (μm)")
    ax.legend(markerscale=1, fontsize=8, loc='best')

    plt.tight_layout()
    out_png = os.path.join(PLOT_FOLDER, f"{slide_id}_allcells.png")
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"Saved plot for {slide_id} → {out_png}")

print("\nAll slides plotted into:", PLOT_FOLDER)
