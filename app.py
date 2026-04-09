import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import zscore

# -----------------------------
# PAGE CONFIG (PREMIUM UI)
# -----------------------------
st.set_page_config(
    page_title="Ubixplorer",
    layout="wide",
)

st.title("🧬 Ubixplorer: Ubiquitin-Based Spatial Explorer for AD–HD Profiling")

st.markdown(
"""
Explore spatial expression patterns, hotspot regions, and biological relevance  
of ubiquitin-related hub genes in Alzheimer's and Huntington's Disease.
"""
)

# -----------------------------
# LOAD DATA
# -----------------------------
expr = pd.read_csv("expression_small.csv")
coords = pd.read_csv("tissue_positions_fixed.csv")

# fix columns
expr = expr.rename(columns={"Unnamed: 0": "BARCODE"})
coords = coords.rename(columns={"barcode": "BARCODE"})

df = coords.merge(expr, on="BARCODE")

# -----------------------------
# SIDEBAR (CONTROL PANEL)
# -----------------------------
st.sidebar.header("⚙️ Controls")

genes = ["MYC","RELA","MYD88","NFKBIA","TLR2","FOXO1"]
gene = st.sidebar.selectbox("Select Gene", genes)

threshold = st.sidebar.slider("Hotspot Threshold (Z-score)", 1.0, 3.0, 2.0)

# -----------------------------
# LAYOUT (2 COLUMNS)
# -----------------------------
col1, col2 = st.columns(2)

# -----------------------------
# SPATIAL EXPRESSION
# -----------------------------
with col1:
    st.subheader(f"📍 Spatial Expression: {gene}")

    fig, ax = plt.subplots(figsize=(5,5))

    sc = ax.scatter(
        df["pxl_col"],
        df["pxl_row"],
        c=df[gene],
        s=5
    )

    ax.invert_yaxis()
    plt.colorbar(sc, ax=ax)

    st.pyplot(fig)

# -----------------------------
# HOTSPOTS
# -----------------------------
df[f"{gene}_z"] = zscore(df[gene])
hotspots = df[df[f"{gene}_z"] > threshold]

with col2:
    st.subheader("🔥 Hotspot Regions")

    fig2, ax2 = plt.subplots(figsize=(5,5))

    ax2.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax2.scatter(hotspots["pxl_col"], hotspots["pxl_row"], color="red", s=6)

    ax2.invert_yaxis()

    st.pyplot(fig2)

# -----------------------------
# INTERPRETATION
# -----------------------------
st.markdown("---")
st.subheader("🧠 Interpretation")

mean_expr = df[gene].mean()
high_expr = hotspots.shape[0]

# expression level
if mean_expr < 0.1:
    level = "low"
elif mean_expr < 0.5:
    level = "moderate"
else:
    level = "high"

# hotspot pattern
if high_expr < 50:
    pattern = "sparse and localized"
elif high_expr < 200:
    pattern = "moderately distributed"
else:
    pattern = "widely distributed"

st.write(f"""
• **Expression Level:** {level}  
• **Hotspot Distribution:** {pattern} ({high_expr} regions)  

This suggests that **{gene}** exhibits {pattern} activity across the tissue,  
indicating its spatial involvement in ubiquitin-mediated processes.
""")

# -----------------------------
# BIOLOGICAL PATHWAYS
# -----------------------------
st.subheader("🔬 Gene-Associated Pathways")

if gene in ["RELA", "NFKBIA"]:
    st.info("NF-κB signaling pathway → central to inflammation and neurodegeneration.")

elif gene in ["MYD88", "TLR2"]:
    st.info("Toll-like receptor / innate immune signaling → indicates immune activation.")

elif gene == "MYC":
    st.info("Cell cycle & transcriptional regulation → linked to cellular stress and proliferation.")

elif gene == "FOXO1":
    st.info("Oxidative stress response & longevity pathways → regulates survival mechanisms.")

# -----------------------------
# FOOTER
# -----------------------------
st.markdown("---")
st.caption("Developed for spatial analysis of ubiquitin-related gene activity in AD vs HD.")
