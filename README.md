# 🧪 Research @ Nationwide Children's Hospital  
**Spatial Analysis of the Tumor Microenvironment in Cancer Immunotherapy**

---

## 📊 Overview

- This research focuses on understanding how the **spatial organization of immune and tumor cells** within the **tumor microenvironment (TME)** affects cancer progression and response to immunotherapy.
- The datasets were derived from **published studies** using **high-resolution spatial single-cell imaging** of tumor tissues from **brain, lung, and head/neck cancers**.
- Developed Python algorithms to:
  - Map cell distributions within tumor tissue slides
  - Identify and **cluster immune cells** (e.g., CD8+ T-cells)
  - Compare spatial patterns between **responders** and **non-responders** to treatment
- Currently co-authoring a manuscript for submission to a peer-reviewed journal.

---

## 🧰 Tools & Technologies

- **Python** (developed in PyCharm)
- Libraries: `pandas`, `matplotlib`, `numpy`, `scipy.io`, `statsmodels`, `scikit-learn`

---

## 🔍 Key Concepts Explained

### 🧫 What Are Immune Cells Doing in Tumors?
- Immune cells like **CD8+ T-cells** are expected to **attack tumor cells**.
- Their **location and clustering** relative to tumor cells can reveal how well the immune system is responding.

### 🧱 What Is the Tumor Microenvironment (TME)?
- The **TME** includes tumor cells, immune cells, blood vessels, and stromal tissue.
- The **spatial layout** of these cells can influence how a tumor grows or responds to treatment.

### 🧪 Responders vs Non-Responders
- Patients are categorized as:
  - **Responders**: their tumor shrinks or stabilizes after **immune checkpoint inhibitor (ICI)** therapy.
  - **Non-Responders**: their tumor does not improve.
- Analyzing these groups helps uncover spatial biomarkers of **treatment success**.

---

## 🧭 Project Pipeline

### 🧬 1. Mapping Out Tumor Slides

We begin by visualizing individual cells from tissue slides using their **x, y spatial coordinates**. These coordinates are sometimes encoded through **lattice correlations** (a method for positioning cells on structured grids).

![Slide Cell Map](BrM_001C2_allcells.png)

📄 [Code for slide plotting](cellMap.py)

---

### 🧪 2. Clustering Immune Cells

CD8+ (cytotoxic T-cells) are clustered using **DBSCAN** to find groups within **30 μm** proximity. Clustering immune cells reveals how they may be "attacking" or surrounding tumor cells.

![CD8 Cluster](tcellclusterView.png)

📄 [Clustering & visualization code](cellMap.py)

---

### 🛡️ 3. Fencing Metric: Measuring Tumor "Surrounding"

To understand immune-tumor interaction, we calculate a **fencing metric**:  
> How well immune cells “fence in” or surround tumor cells.

Comparing **responders** and **non-responders** with this metric helps assess whether **immune cell positioning** relates to successful treatment.

![Alt Mac Hist](Alt_MAC_hist.png)

📄 [Statistical analysis & fencing calculations](NRvR.py)

---

## 📏 Fencing Metric Parameters

| Parameter              | Description                                                  |
|------------------------|--------------------------------------------------------------|
| `N_SHUFFLES = 5000`    | Number of random permutations for statistical significance    |
| `DBSCAN_EPS = 30`      | Max distance between clustered immune cells (μm)             |
| `DBSCAN_MIN_SAMPLES=3` | Minimum number of immune cells to form a valid cluster       |
| `MIN_CELL_COUNT = 20`  | Minimum cell count required on a slide for inclusion         |
| `ALPHA = 0.05`         | Significance level for hypothesis testing                    |

These constraints ensure **statistical rigor** and eliminate noise from low-cell-count slides.

---

