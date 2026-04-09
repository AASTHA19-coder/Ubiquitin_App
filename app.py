# -----------------------------
# TABS (MAIN UI UPGRADE)
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
# OVERVIEW TAB
# -----------------------------
with tabs[0]:
    st.subheader("🧬 Study Overview")

    st.info("""
This tool is developed as part of our study integrating:

• Differential gene expression analysis (AD vs HD vs Control)  
• Identification of ubiquitin-related hub genes  
• Spatial transcriptomics (Visium)  

We map where these genes are active in brain tissue  
to understand localized neurodegenerative processes.
""")

# -----------------------------
# EXPRESSION TAB
# -----------------------------
with tabs[1]:
    st.subheader(f"📍 Spatial Expression: {gene}")

    fig, ax = plt.subplots()
    sc = ax.scatter(df["pxl_col"], df["pxl_row"], c=df[gene], s=5)

    ax.invert_yaxis()
    plt.colorbar(sc, ax=ax)

    st.pyplot(fig)

    st.caption("Each point = spatial spot in brain tissue")

# -----------------------------
# HOTSPOTS TAB
# -----------------------------
with tabs[2]:
    st.subheader("🔥 Hotspot Regions")

    df[f"{gene}_z"] = zscore(df[gene])
    hotspots = df[df[f"{gene}_z"] > threshold]

    fig, ax = plt.subplots()
    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(hotspots["pxl_col"], hotspots["pxl_row"], color="red", s=6)

    ax.invert_yaxis()
    st.pyplot(fig)

    st.write("Red regions = high expression zones (Z-score based)")

# -----------------------------
# RISK TAB
# -----------------------------
with tabs[3]:
    st.subheader("🚨 Top Risk Regions")

    top_regions = df.sort_values(f"{gene}_z", ascending=False).head(50)

    fig, ax = plt.subplots()
    ax.scatter(df["pxl_col"], df["pxl_row"], color="lightgrey", s=2)
    ax.scatter(top_regions["pxl_col"], top_regions["pxl_row"], color="yellow", s=10)

    ax.invert_yaxis()
    st.pyplot(fig)

    st.write("Top regions = highest gene activity")

# -----------------------------
# INTERPRETATION TAB
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

    st.write(f"""
In our study, **{gene}** shows **{level} expression**.

We identified **{count} high-activity regions**,  
indicating localized involvement in neurodegeneration.

These findings align with our DEG-based hub gene analysis.
""")

# -----------------------------
# PATHWAYS TAB
# -----------------------------
with tabs[5]:
    st.subheader("🔬 Gene Pathways")

    if gene in ["RELA", "NFKBIA"]:
        st.success("NF-κB pathway → inflammation")

    elif gene in ["MYD88", "TLR2"]:
        st.success("Innate immune signaling")

    elif gene == "MYC":
        st.success("Cell cycle regulation")

    elif gene == "FOXO1":
        st.success("Oxidative stress response")
