import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------
# PAGE CONFIG
# -------------------------------
st.set_page_config(layout="wide")

# -------------------------------
# LOAD DATA
# -------------------------------
expr = pd.read_csv("expression_small.csv")
coords = pd.read_csv("tissue_positions_fixed.csv", header=None)

# Fix columns (Visium format)
coords.columns = ["BARCODE", "IN_TISSUE", "ARRAY_ROW", "ARRAY_COL", "PXL_ROW", "PXL_COL"]

# Rename expression barcode column
expr = expr.rename(columns={"Unnamed: 0": "BARCODE"})

# Merge
df = coords.merge(expr, on="BARCODE")

# -------------------------------
# CUSTOM STYLING (Premium UI)
# -------------------------------
st.markdown("""
<style>
body {
    font-family: 'Segoe UI', sans-serif;
}
.big-title {
    font-size:40px;
    font-weight:700;
    margin-bottom:0px;
}
.subtitle {
    font-size:18px;
    color:gray;
    margin-bottom:30px;
}
.section-title {
    font-size:24px;
    font-weight:600;
    margin-top:30px;
}
</style>
""", unsafe_allow_html=True)

# -------------------------------
# HEADER WITH LOGO
# -------------------------------
col_logo, col_title = st.columns([1, 4])

with col_logo:
    st.image("images/UbiXplorer logo desi.png", width=160)

with col_title:
    st.markdown('<div class="big-title">UbiXplorer</div>', unsafe_allow_html=True)
    st.markdown('<div class="subtitle">Ubiquitin-based spatial profiling for Alzheimer’s vs Huntington’s disease</div>', unsafe_allow_html=True)

# -------------------------------
# SIDEBAR CONTROLS
# -------------------------------
st.sidebar.header("Controls")

gene = st.sidebar.selectbox(
    "Select Gene",
    ["MYC", "RELA", "MYD88", "NFKBIA", "TLR2", "FOXO1"]
)

threshold = st.sidebar.slider("Hotspot Threshold (Z-score)", 1.0, 3.0, 2.0)

# -------------------------------
# SPATIAL EXPRESSION PLOT
# -------------------------------
st.markdown('<div class="section-title">Spatial Expression</div>', unsafe_allow_html=True)

fig, ax = plt.subplots()

sc = ax.scatter(
    df["PXL_COL"],
    df["PXL_ROW"],
    c=df[gene],
    s=10
)

plt.colorbar(sc)
ax.invert_yaxis()

st.pyplot(fig)

# -------------------------------
# HOTSPOT DETECTION
# -------------------------------
st.markdown('<div class="section-title">Hotspot Regions</div>', unsafe_allow_html=True)

z = (df[gene] - df[gene].mean()) / df[gene].std()
hotspots = z > threshold

fig2, ax2 = plt.subplots()

ax2.scatter(
    df["PXL_COL"],
    df["PXL_ROW"],
    c=hotspots.astype(int),
    s=10
)

ax2.invert_yaxis()

st.pyplot(fig2)

# -------------------------------
# INTERPRETATION
# -------------------------------
st.markdown('<div class="section-title">Interpretation</div>', unsafe_allow_html=True)

num_hotspots = hotspots.sum()

st.write(f"• Total hotspot regions: {num_hotspots}")

if num_hotspots < 50:
    st.write("• Pattern: Sparse expression")
else:
    st.write("• Pattern: Widespread activation")

# Gene-specific biology
if gene in ["RELA", "NFKBIA"]:
    st.info(" NF-κB pathway activation → inflammation-related signaling")

elif gene == "MYC":
    st.info(" MYC → transcriptional regulation and neuronal stress")

elif gene == "MYD88":
    st.info(" MYD88 → innate immune signaling")

elif gene == "TLR2":
    st.info(" TLR2 → Toll-like receptor pathway (immune activation)")

elif gene == "FOXO1":
    st.info(" FOXO1 → oxidative stress and apoptosis regulation")

# -------------------------------
# FOOTER
# -------------------------------
st.markdown("---")
st.caption("Human brain Visium dataset • Spatial + quantitative analysis of ubiquitin-related hub genes")
