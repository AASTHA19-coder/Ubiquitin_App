import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

st.title("🧬 Ubiquitin Spatial Disease Explorer")

# -----------------------------
# LOAD DATA
# -----------------------------
expr = pd.read_csv("expression_small.csv")
coords = pd.read_csv("tissue_positions.csv")

# -----------------------------
# FIX COLUMN NAMES (ROBUST)
# -----------------------------
# uppercase everything
expr.columns = expr.columns.str.upper()
coords.columns = coords.columns.str.upper()

# fix barcode column
if "UNNAMED: 0" in expr.columns:
    expr = expr.rename(columns={"UNNAMED: 0": "BARCODE"})

if "BARCODE" not in coords.columns:
    coords = coords.rename(columns={coords.columns[0]: "BARCODE"})

# -----------------------------
# MERGE
# -----------------------------
df = coords.merge(expr, on="BARCODE")

# -----------------------------
# AUTO-DETECT COORDINATES
# -----------------------------
col_candidates = [c for c in df.columns if "COL" in c.upper()]
row_candidates = [c for c in df.columns if "ROW" in c.upper()]

pxl_col = col_candidates[0]
pxl_row = row_candidates[0]

# -----------------------------
# SELECT GENE
# -----------------------------
genes = ["MYC","RELA","MYD88","NFKBIA","TLR2","FOXO1"]
gene = st.selectbox("Select Gene", genes)

# -----------------------------
# SPATIAL EXPRESSION
# -----------------------------
st.subheader(f"{gene} Spatial Expression")

fig, ax = plt.subplots(figsize=(5,5))

sc = ax.scatter(
    df[pxl_col],
    df[pxl_row],
    c=df[gene],
    s=5
)

ax.invert_yaxis()
plt.colorbar(sc, ax=ax)

st.pyplot(fig)

# -----------------------------
# HOTSPOT DETECTION
# -----------------------------
st.subheader("Hotspot Detection")

threshold = st.slider("Z-score Threshold", 1.0, 3.0, 2.0)

df[f"{gene}_z"] = zscore(df[gene])
hotspots = df[df[f"{gene}_z"] > threshold]

fig2, ax2 = plt.subplots(figsize=(5,5))

ax2.scatter(df[pxl_col], df[pxl_row], color="lightgrey", s=2)
ax2.scatter(hotspots[pxl_col], hotspots[pxl_row], color="red", s=6)

ax2.invert_yaxis()

st.pyplot(fig2)

# -----------------------------
# ABOUT
# -----------------------------
st.subheader("About")
st.write(
    "Interactive visualization of ubiquitin-related hub gene expression and spatial hotspots "
    "in brain tissue datasets."
)
