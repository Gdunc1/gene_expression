import streamlit as st

st.title("Gene Expression Analysis and Data Exploration App")
st.subheader("Data Visualization Final Project - FALL 2023")
st.subheader("Developed by: Gerard Duncan")
st.subheader("Data Source: Broad Institute depmap")

st.write("""
    Welcome to my application designed for gene expression analysis and target identification using the Cancer Cell Line Encyclopedia (CCLE) dataset, which includes over 1000 cell line models. This app facilitates the analysis of differential gene expression (DGE) and provides an interactive platform for exploring genetic mutation data.
    
    #### Features
    - **Data Exploration**: Explore mutation data across different diseases. Visualize the most common mutated genes and analyze disease-specific trends.
    - **DGE Analysis**: Perform differential gene expression analysis, focusing on NSCLC cell lines. Users can select specific genes and visualize the analysis results.
    - **Interactive Volcano Plot**: Upload your DGE analysis results and visualize them in an interactive volcano plot, or explore demo results for insights.

    Navigate through the app using the sidebar to access different functionalities. Start by selecting a page from the options provided.
    """)
