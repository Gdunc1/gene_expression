import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

def create_volcano_plot(dataframe):
    df = dataframe.copy()
    df['gene'] = df.index
    df["-log10(pvalue)"] = -np.log10(df.pvalue)
    df["-log10(padj)"] = -np.log10(df.padj)

    fig = px.scatter(
        df, width=800, height=800,
        x="log2FoldChange",
        y="-log10(padj)",
        color="-log10(pvalue)",
        hover_data=["gene"],
        labels={"log2FoldChange": "log2 Fold Change", "-log10(padj)": "-log10(p adj)"},
        title="Interactive Volcano Plot",
        color_continuous_scale=["blue", "red"]  # Custom color scale: blue to red
    )
    fig.update_traces(marker=dict(size=8, opacity=0.5))
    return fig

# Load the demo file
demo_file = 'data/DGE_results_demo.csv'
demo_data = pd.read_csv(demo_file, index_col=0)

# Streamlit layout for volcano plot page
st.title("Interactive Volcano Plot of DGE Results")
st.write("""
This page provides an interactive visualization tool - the Volcano Plot - to explore Differential Gene Expression (DGE) results. 

## What is a Volcano Plot?
A volcano plot is a type of scatter plot used to quickly identify changes in large datasets composed of replicate data. It plots significance versus fold-change on the y and x axes, respectively.

## How to Use:
1. **Upload Your DGE Results**: Use the uploader to input your DGE analysis results in CSV format. The file should contain columns for p-values, adjusted p-values (padj), and log2 fold change.
2. **Interact with the Plot**: Once uploaded, the volcano plot will be displayed. You can hover over points to see gene names and their respective values.
3. **Explore and Analyze**: Identify genes with significant differential expression. Genes with high fold changes and low p-values will be highlighted.

Bellow is an example of the analysis output showing Volcano Plot of DGE results for TP53-mut vs TP53-wt NSCLC cell lines.
""")

uploaded_file = st.sidebar.file_uploader("Upload DGE Results CSV", type=["csv"])

# Use uploaded data if available, else use demo data
if uploaded_file is not None:
    dge_results = pd.read_csv(uploaded_file, index_col=0)
else:
    dge_results = demo_data

volcano_plot = create_volcano_plot(dge_results)
st.plotly_chart(volcano_plot)
