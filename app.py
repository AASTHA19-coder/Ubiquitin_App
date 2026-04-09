import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

# -----------------------------
# PAGE CONFIG
# -----------------------------
st.set_page_config(page_title="Ubixplorer", layout="wide")

# -----------------------------
# CUSTOM CSS (FONTS + STYLE)
# -----------------------------
st.markdown("""
<style>
body {
    font-family: 'Segoe UI', sans-serif;
}
h1 {
    font-size: 50px !important;
    font-weight: 700;
}
h2 {
    font-size: 30px !important;
    font-weight: 600;
}
h3 {
    font-size: 25px !important;
    font-weight: 600;
}
p {
    font-size: 18px;
}
.block-container {
    padding-top: 2rem;
}
</style>
""", unsafe_allow_html=True)

# -----------------------------
# TITLE
# -----------------------------
st.markdown("## 🧬 UbiXplorer")
st.markdown("### Ubiquitin-Based Spatial Explorer for AD–HD Profiling")

st.markdown("""
This platform integrates differential expression analysis, hub gene identification, and spatial transcriptomics to investigate ubiquitin-related mechanisms in neurodegenerative diseases.
""")

# -----------------------------
# LOAD DATA
# -----------------------------
expr = pd.read_csv("expression_small.csv")
coords = pd.read_csv("tissue_positions_fixed.csv")

expr = expr.rename(columns={"Unnamed: 0": "BARCODE"})
coords.columns = ["BARCODE", "IN_TISSUE", "ARRAY_ROW", "ARRAY_COL", "pxl_row", "pxl_col"]

df = coords.merge(expr, on="BARCODE")

# -----------------------------
# SIDEBAR
# -----------------------------
st.sidebar.markdown("## ⚙️ Controls")

genes = ["MYC","RELA","MYD88","NFKBIA","TLR2","FOXO1"]
gene = st.sidebar.selectbox("Select Gene", genes)

threshold = st.sidebar.slider("Hotspot Threshold (Z-score)", 1.0, 3.0, 2.0)

# -----------------------------
# COMPUTATION
# -----------------------------
df[f"{gene}_z"] = zscore(df[gene].fillna(0))

hotspots = df[df[f"{gene}_z"] > threshold]
top_regions = df.sort_values(f"{gene}_z", ascending=False).head(50)

# -----------------------------
# TABS
# -----------------------------
tabs = st.tabs([
    "Overview",
    "Expression",
    "Hotspots",
    "Risk Regions",
    "Interpretation",
    "Pathways"
])

# -----------------------------
# OVERVIEW
# -----------------------------
with tabs[0]:
    st.markdown("### **Study Overview**")

    st.markdown("""
- Comparative analysis performed across:
  - AD vs Control  
  - AD vs HD  
  - HD vs Control  

- Identification of **ubiquitin-related hub genes**

- Projection onto **human brain Visium spatial transcriptomics data**

- Enables **spatial localization of disease-associated molecular activity**
""")

# -----------------------------
# EXPRESSION
# -----------------------------
with tabs[1]:
    st.markdown(f"### **Spatial Expression: {gene}**")

    fig, ax = plt.subplots()

    sc = ax.scatter(df["pxl_col"], df["pxl_row"], c=df[gene], s=6)

    #ax.set_xlabel("Tissue X Coordinate")
    #ax.set_ylabel("Tissue Y Coordinate")
    #ax.invert_yaxis()

    plt.colorbar(sc, ax=ax)

    st.pyplot(fig)

    st.markdown("""
**Interpretation:**  
x-axis and y-axis represnt the tissue coordinates.

Each point represents a spatial location on the human brain tissue section.  

Color intensity reflects gene expression levels.
""")

# -----------------------------
# HOTSPOTS
# -----------------------------
with tabs[2]:
    st.markdown("### **Hotspot Regions (Z-score based)**")

    fig, ax = plt.subplots()

    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(hotspots["pxl_col"], hotspots["pxl_row"], color="red", s=8)

    #ax.set_xlabel("Tissue X Coordinate")
    #ax.set_ylabel("Tissue Y Coordinate")
    #ax.invert_yaxis()

    st.pyplot(fig)

    st.markdown(f"**Detected Hotspots:** {len(hotspots)} regions"
               "x-axis and y-axis represnt the tissue coordinates.")

# -----------------------------
# RISK REGIONS
# -----------------------------
with tabs[3]:
    st.markdown("### **Top Risk Regions**")

    fig, ax = plt.subplots()

    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(top_regions["pxl_col"], top_regions["pxl_row"], color="yellow", s=12)

    #ax.set_xlabel("Tissue X Coordinate")
    #ax.set_ylabel("Tissue Y Coordinate")
    #ax.invert_yaxis()

    st.pyplot(fig)

    st.markdown("""
x-axis and y-axis represnt the tissue coordinates.

Top 50 regions with highest expression levels → potential disease-relevant zones.
""")

# -----------------------------
# INTERPRETATION
# -----------------------------
with tabs[4]:
    st.markdown("### **Biological Interpretation**")

    st.markdown(f"""
The gene **{gene}** was identified as a hub gene in AD–HD comparative analysis.  
Its spatial distribution highlights localized functional activity.
""")

    st.markdown("#### **Disease Context**")

    if gene == "MYC":
        st.write("Associated with **AD vs Control** → cell cycle dysregulation")

    elif gene in ["RELA", "NFKBIA"]:
        st.write("Shared in AD & HD → **NF-κB inflammatory signaling**")

    elif gene in ["MYD88", "TLR2"]:
        st.write("Innate immune activation → **neuroinflammation**")

    elif gene == "FOXO1":
        st.write("Oxidative stress → **neuronal survival mechanisms**")

# -----------------------------
# PATHWAYS
# -----------------------------
with tabs[5]:
    st.markdown("### **Gene-Associated Pathways**")

    if gene in ["RELA", "NFKBIA"]:
        st.success("NF-κB signaling pathway → inflammation")

    elif gene in ["MYD88", "TLR2"]:
        st.success("Toll-like receptor signaling → immune response")

    elif gene == "MYC":
        st.success("Cell cycle regulation")

    elif gene == "FOXO1":
        st.success("Oxidative stress & apoptosis")

# -----------------------------
# FOOTER
# -----------------------------
st.markdown("---")
st.caption("Ubixplorer | Spatial mapping of ubiquitin-related hub genes in AD and HD")
