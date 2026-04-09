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
This interactive platform is developed from our study integrating differential expression analysis,  
hub gene identification, and spatial transcriptomics to explore ubiquitin-related mechanisms  
in neurodegenerative diseases (Alzheimer’s Disease and Huntington’s Disease).
""")

# -----------------------------
# LOAD DATA
# -----------------------------
expr = pd.read_csv("expression_small.csv")
coords = pd.read_csv("tissue_positions_fixed.csv")

# Fix column names
expr = expr.rename(columns={"Unnamed: 0": "BARCODE"})

coords.columns = ["BARCODE", "IN_TISSUE", "ARRAY_ROW", "ARRAY_COL", "pxl_row", "pxl_col"]

# Merge
df = coords.merge(expr, on="BARCODE")

# -----------------------------
# SIDEBAR
# -----------------------------
st.sidebar.header("⚙️ Controls")

genes = ["MYC","RELA","MYD88","NFKBIA","TLR2","FOXO1"]
gene = st.sidebar.selectbox("Select Gene", genes)

threshold = st.sidebar.slider("Hotspot Threshold (Z-score)", 1.0, 3.0, 2.0)

# -----------------------------
# COMPUTATIONS
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
    col1, col2 = st.columns([1,5])

    with col1:
        st.image("images/overview.png", width=80)

    with col2:
        st.subheader("Study Overview")

        st.info("""
We performed comparative analysis across:

• AD vs Control  
• AD vs HD  
• HD vs Control  

From these, we identified **ubiquitin-related hub genes**.

These hub genes were then projected onto **human brain Visium spatial transcriptomics data**  
to understand their **spatial localization and activity patterns**.

This platform enables exploration of how these disease-associated genes behave in real tissue.
""")

# -----------------------------
# EXPRESSION
# -----------------------------
with tabs[1]:
    col1, col2 = st.columns([1,5])

    with col1:
        st.image("images/expression.png", width=80)

    with col2:
        st.subheader(f"Spatial Expression: {gene}")

    fig, ax = plt.subplots()

    sc = ax.scatter(df["pxl_col"], df["pxl_row"], c=df[gene], s=6)

    ax.set_xlabel("Tissue X Coordinate")
    ax.set_ylabel("Tissue Y Coordinate")
    ax.invert_yaxis()

    plt.colorbar(sc, ax=ax)

    st.pyplot(fig)

    st.caption("""
Each point represents a spatial spot on a human brain tissue section (Visium platform).  
Color intensity indicates expression level of the selected gene.
""")

# -----------------------------
# HOTSPOTS
# -----------------------------
with tabs[2]:
    col1, col2 = st.columns([1,5])

    with col1:
        st.image("images/hotspot.png", width=80)

    with col2:
        st.subheader("Hotspot Regions (Z-score based)")

    fig, ax = plt.subplots()

    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(hotspots["pxl_col"], hotspots["pxl_row"], color="red", s=8)

    ax.set_xlabel("Tissue X Coordinate")
    ax.set_ylabel("Tissue Y Coordinate")
    ax.invert_yaxis()

    st.pyplot(fig)

    st.write(f"Total hotspot regions detected: {len(hotspots)}")

# -----------------------------
# RISK REGIONS
# -----------------------------
with tabs[3]:
    col1, col2 = st.columns([1,5])

    with col1:
        st.image("images/risk.png", width=80)

    with col2:
        st.subheader("Top Risk Regions")

    fig, ax = plt.subplots()

    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(top_regions["pxl_col"], top_regions["pxl_row"], color="yellow", s=12)

    ax.set_xlabel("Tissue X Coordinate")
    ax.set_ylabel("Tissue Y Coordinate")
    ax.invert_yaxis()

    st.pyplot(fig)

    st.write("Top 50 regions with highest expression (potential disease-relevant zones).")

# -----------------------------
# INTERPRETATION
# -----------------------------
with tabs[4]:
    col1, col2 = st.columns([1,5])

    with col1:
        st.image("images/interpretation.png", width=80)

    with col2:
        st.subheader("Biological Interpretation")

    st.write(f"""
The gene **{gene}** was identified as a hub gene in our differential analysis.

Its spatial pattern suggests localized biological activity within the brain tissue.
""")

    st.markdown("### Disease Context")

    if gene == "MYC":
        st.write("Strongly associated with **AD vs Control** → cell cycle dysregulation.")

    elif gene in ["RELA", "NFKBIA"]:
        st.write("Shared in AD and HD → **NF-κB inflammatory pathway activation**.")

    elif gene in ["MYD88", "TLR2"]:
        st.write("Innate immune signaling → **neuroinflammation in both diseases**.")

    elif gene == "FOXO1":
        st.write("Stress response → **neuronal survival and degeneration balance**.")

# -----------------------------
# PATHWAYS
# -----------------------------
with tabs[5]:
    col1, col2 = st.columns([1,5])

    with col1:
        st.image("images/pathways.png", width=80)

    with col2:
        st.subheader("Gene-Associated Pathways")

    if gene in ["RELA", "NFKBIA"]:
        st.success("NF-κB signaling → inflammation")

    elif gene in ["MYD88", "TLR2"]:
        st.success("Toll-like receptor / innate immune signaling")

    elif gene == "MYC":
        st.success("Cell cycle regulation and proliferation")

    elif gene == "FOXO1":
        st.success("Oxidative stress response and apoptosis")

# -----------------------------
# FOOTER
# -----------------------------
st.markdown("---")
st.caption("Ubixplorer | Spatial mapping of ubiquitin-related hub genes in AD and HD")
