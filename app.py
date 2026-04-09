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
Explore spatial expression patterns, hotspot regions, and biological relevance  
of ubiquitin-related hub genes identified from AD vs HD analysis.
""")

# -----------------------------
# ABOUT SECTION
# -----------------------------
with st.expander("📖 About this Analysis"):
    st.markdown("""
### 🧬 Dataset
Human brain spatial transcriptomics data (10x Genomics Visium)

### 📍 Axes
- X-axis → spatial column (pxl_col)  
- Y-axis → spatial row (pxl_row)  

Each point represents a spatial spot on the tissue.

### 🔬 Workflow
1. DEG analysis (AD vs Control, HD vs Control, AD vs HD)  
2. Identification of ubiquitin-related hub genes  
3. Spatial mapping on brain tissue  
4. Hotspot detection using Z-score  

### 🎯 Purpose
To identify spatial activation patterns of disease-associated genes  
and highlight potential regions involved in neurodegeneration.
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

# -----------------------------
# PLOTS
# -----------------------------
col1, col2 = st.columns(2)

# Spatial Expression
with col1:
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

# Hotspots
df[f"{gene}_z"] = zscore(df[gene])
hotspots = df[df[f"{gene}_z"] > threshold]

with col2:
    st.subheader("🔥 Hotspot Regions")

    fig2, ax2 = plt.subplots()

    ax2.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax2.scatter(hotspots["pxl_col"], hotspots["pxl_row"], color="red", s=6)

    ax2.invert_yaxis()

    st.pyplot(fig2)

st.caption("📍 X and Y axes represent spatial coordinates on the brain tissue (Visium platform).")

# -----------------------------
# TOP RISK REGIONS
# -----------------------------
st.markdown("---")
st.subheader("🚨 Top Risk Regions")

top_regions = df.sort_values(f"{gene}_z", ascending=False).head(50)

fig3, ax3 = plt.subplots()

ax3.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
ax3.scatter(top_regions["pxl_col"], top_regions["pxl_row"], color="yellow", s=10)

ax3.invert_yaxis()

st.pyplot(fig3)

st.write("Top regions represent highest gene activity → potential disease-relevant zones.")

# -----------------------------
# INTERPRETATION
# -----------------------------
st.markdown("---")
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
    pattern = "localized"
elif count < 200:
    pattern = "moderately distributed"
else:
    pattern = "widely distributed"

st.write(f"""
• Expression level: **{level}**  
• Spatial pattern: **{pattern}**  
• Top active regions: **{count}**

This suggests that **{gene}** exhibits {pattern} activity across the tissue,  
indicating its involvement in ubiquitin-mediated neurodegenerative processes.
""")

# -----------------------------
# DISEASE RELEVANCE (KEY ADDITION)
# -----------------------------
st.subheader("🧬 Disease Relevance (From DEG Analysis)")

if gene == "MYC":
    st.write("🔴 Strongly associated with **Alzheimer’s Disease (AD vs Control)**")

elif gene in ["RELA", "NFKBIA", "MYD88", "TLR2"]:
    st.write("🟠 Shared hub genes across **AD and HD**, linked to inflammation and ubiquitin pathways")

elif gene == "FOXO1":
    st.write("🟡 Associated with stress response in neurodegeneration")

else:
    st.write("General ubiquitin-related regulatory gene")

# -----------------------------
# PATHWAYS
# -----------------------------
st.subheader("🔬 Gene-Associated Pathways")

if gene in ["RELA", "NFKBIA"]:
    st.info("NF-κB signaling pathway → inflammation and neurodegeneration")

elif gene in ["MYD88", "TLR2"]:
    st.info("Toll-like receptor pathway → innate immune activation")

elif gene == "MYC":
    st.info("Cell cycle regulation → cellular stress and transcription control")

elif gene == "FOXO1":
    st.info("Oxidative stress response → longevity and survival pathways")

# -----------------------------
# FOOTER
# -----------------------------
st.markdown("---")
st.caption("Ubixplorer | Spatial analysis of ubiquitin-related hub genes in AD and HD")
