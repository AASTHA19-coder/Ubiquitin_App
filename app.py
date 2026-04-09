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

# merge
df = coords.merge(expr, on="BARCODE")

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
    df["PXL_COL"],
    df["PXL_ROW"],
    c=df[gene],
    s=5
)

ax.invert_yaxis()
plt.colorbar(sc, ax=ax)

st.pyplot(fig)

# -----------------------------
# HOTSPOTS
# -----------------------------
st.subheader("Hotspot Detection")

threshold = st.slider("Z-score Threshold", 1.0, 3.0, 2.0)

df[f"{gene}_z"] = zscore(df[gene])
hotspots = df[df[f"{gene}_z"] > threshold]

fig2, ax2 = plt.subplots(figsize=(5,5))

ax2.scatter(df["PXL_COL"], df["PXL_ROW"], color="lightgrey", s=2)
ax2.scatter(hotspots["PXL_COL"], hotspots["PXL_ROW"], color="red", s=6)

ax2.invert_yaxis()

st.pyplot(fig2)

# -----------------------------
# ABOUT
# -----------------------------
st.subheader("About")
st.write("Interactive exploration of ubiquitin-related gene expression in spatial transcriptomics data.")