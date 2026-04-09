import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

# -----------------------------
# PAGE CONFIG
# -----------------------------
st.set_page_config(page_title="Ubixplorer", layout="wide")

st.title("🧬 Ubixplorer: Ubiquitin-Based Spatial Explorer for AD–HD Profiling")

st.markdown("""
An interactive platform developed from our study integrating DEG analysis,  
hub gene identification, and spatial transcriptomics to understand  
ubiquitin-related mechanisms in neurodegeneration.
""")

# -----------------------------
# LOAD DATA
# -----------------------------
expr = pd.read_csv("expression_small.csv")
coords = pd.read_csv("tissue_positions_fixed.csv")

expr = expr.rename(columns={"Unnamed: 0": "BARCODE"})
coords = coords.rename(columns={"barcode": "BARCODE"})

df = coords.merge(expr, on="BARCODE")

# -----------------------------
# SIDEBAR
# -----------------------------
st.sidebar.header("⚙️ Controls")

genes = ["MYC","RELA","MYD88","NFKBIA","TLR2","FOXO1"]
gene = st.sidebar.selectbox("Select Gene", genes)

threshold = st.sidebar.slider("Hotspot Threshold (Z-score)", 1.0, 3.0, 2.0)

# Precompute
df[f"{gene}_z"] = zscore(df[gene])
hotspots = df[df[f"{gene}_z"] > threshold]
top_regions = df.sort_values(f"{gene}_z", ascending=False).head(50)

# -----------------------------
# TABS UI
# -----------------------------
tabs = st.tabs([
    "🧬 Overview",
    "📍 Expression",
    "🔥 Hotspots",
    "🚨 Risk Regions",
    "🧠 Interpretation",
    "🔬 Pathways"
])

# -----------------------------
# OVERVIEW
# -----------------------------
with tabs[0]:
    st.subheader("🧬 Study Overview")

    st.info("""
This application is derived from our analysis of Alzheimer's Disease (AD)  
and Huntington’s Disease (HD) using transcriptomic data.

We first performed differential expression and network analysis to identify  
ubiquitin-related hub genes. These genes were then mapped onto spatial  
transcriptomics data (Visium) to examine their localization in brain tissue.

This enables exploration of how disease-associated genes are spatially organized.
""")

# -----------------------------
# EXPRESSION
# -----------------------------
with tabs[1]:
    st.subheader(f"📍 Spatial Expression: {gene}")

    fig, ax = plt.subplots()

    sc = ax.scatter(
        df["pxl_col"],
        df["pxl_row"],
        c=df[gene],
        s=5
    )

    ax.invert_yaxis()
    plt.colorbar(sc, ax=ax)

    st.pyplot(fig)

    st.caption("X and Y axes represent spatial coordinates on the Visium brain tissue slide.")

# -----------------------------
# HOTSPOTS
# -----------------------------
with tabs[2]:
    st.subheader("🔥 Hotspot Regions")

    fig, ax = plt.subplots()

    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(hotspots["pxl_col"], hotspots["pxl_row"], color="red", s=6)

    ax.invert_yaxis()
    st.pyplot(fig)

    st.write("""
Hotspots are defined using Z-score normalization.  
Red regions indicate significantly elevated gene expression.
""")

# -----------------------------
# RISK REGIONS
# -----------------------------
with tabs[3]:
    st.subheader("🚨 Top Risk Regions")

    fig, ax = plt.subplots()

    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(top_regions["pxl_col"], top_regions["pxl_row"], color="yellow", s=10)

    ax.invert_yaxis()
    st.pyplot(fig)

    st.write("""
Top-ranked regions represent highest expression levels  
and may correspond to disease-relevant zones.
""")

# -----------------------------
# INTERPRETATION
# -----------------------------
with tabs[4]:
    st.subheader("🧠 Interpretation")

    mean_expr = df[gene].mean()
    count = len(top_regions)

    if mean_expr < 0.1:
        level = "low"
    elif mean_expr < 0.5:
        level = "moderate"
    else:
        level = "high"

    if count < 50:
        pattern = "highly localized"
    elif count < 200:
        pattern = "regionally distributed"
    else:
        pattern = "widely distributed"

    st.write(f"""
In our analysis, **{gene}** demonstrates **{level} expression** across the tissue.

We observe a **{pattern} spatial pattern**, with {count} high-activity regions.

This is consistent with our DEG and network-based identification of  
**{gene}** as a ubiquitin-related hub gene in neurodegenerative conditions.
""")

    # Disease relevance
    st.markdown("### 🧬 Disease Context")

    if gene == "MYC":
        st.write("Identified as a key regulator in **AD vs Control**, suggesting involvement in Alzheimer's pathology.")

    elif gene in ["RELA", "NFKBIA", "MYD88", "TLR2"]:
        st.write("These genes were consistently observed across **AD and HD**, indicating shared inflammatory mechanisms.")

    elif gene == "FOXO1":
        st.write("Associated with stress response and neuronal survival pathways.")

# -----------------------------
# PATHWAYS
# -----------------------------
with tabs[5]:
    st.subheader("🔬 Gene-Associated Pathways")

    if gene in ["RELA", "NFKBIA"]:
        st.success("NF-κB signaling → central to inflammation and neurodegeneration")

    elif gene in ["MYD88", "TLR2"]:
        st.success("Toll-like receptor pathway → innate immune activation")

    elif gene == "MYC":
        st.success("Cell cycle and transcriptional regulation")

    elif gene == "FOXO1":
        st.success("Oxidative stress and longevity pathways")

# -----------------------------
# FOOTER
# -----------------------------
st.markdown("---")
st.caption("Ubixplorer | Spatial mapping of ubiquitin-related hub genes in AD and HD")
