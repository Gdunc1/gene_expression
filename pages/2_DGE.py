import streamlit as st
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from io import StringIO
import base64

NSCLC_celllines = pd.read_csv('data/NSCLC_celllines_filter.csv')
mutation = pd.read_csv('data/OmicsSomaticMutations_filter.csv', index_col=0)
mutation_selection = mutation.HugoSymbol.unique()

demo_file = 'data/DGE_results_demo.csv'
demo_data = pd.read_csv(demo_file, index_col=0)

# Function to create a download link
def download_link(object_to_download, download_filename, download_link_text):
    """
    Generates a link to download the given object_to_download.
    """
    if isinstance(object_to_download, pd.DataFrame):
        object_to_download = object_to_download.to_csv(index=False)

    buffer = StringIO()
    buffer.write(object_to_download)
    buffer.seek(0)
    return st.markdown(f"[{download_link_text}](data:text/csv;charset=utf-8,{buffer.getvalue()})", unsafe_allow_html=True)

def get_table_download_link(df, filename, text):
    csv = df.to_csv(index=True)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}">{text}</a>'
    return href

# DGE Analysis Function
def perform_dge_analysis():
    rnaseq = pd.read_csv('data/RNA_Non_Small_Cell_Lung_Cancer.csv', index_col=0)
    
    mutation_filter = mutation[mutation.HugoSymbol == 'TP53']
    NSCLC_celllines['TP53_status'] = NSCLC_celllines.ModelID.isin(mutation_filter.ModelID)

    # rename columns in rnaseq_filter keeping only the first part of the string before the first space
    rnaseq.columns = rnaseq.columns.str.split(' ').str[0]
    counts_df = rnaseq
    # filter counts_df keeping only ModelID that are in NSCLC_celllines
    counts_df = counts_df[counts_df.index.isin(NSCLC_celllines.ProfileID)]

    # make sure that count_df is a mtrix that only contains integers
    counts_df = counts_df.astype(int)

    # create metadata dataframe with ModelID as an index and TP53_status
    metadata = NSCLC_celllines[['ProfileID', 'TP53_status']].set_index('ProfileID')

    # replace True with 1 and False with 0 in metadata
    metadata.TP53_status = metadata.TP53_status.replace({True: 'TP53_mut', False: 'TP53_wt'})

    # filter index of metadata keeping only ModelID that are in counts_df
    metadata = metadata[metadata.index.isin(counts_df.index)]

    # Order metadata by counts_df index
    metadata = metadata.reindex(counts_df.index)

    samples_to_keep = ~metadata.TP53_status.isna()
    counts_df = counts_df.loc[samples_to_keep]
    metadata = metadata.loc[samples_to_keep]

    # keep the first column among duplicated columns
    counts_df = counts_df.loc[:, ~counts_df.columns.duplicated()]

    dds = DeseqDataSet(
        counts=counts_df,
        metadata=metadata,
        design_factors="TP53_status",
        refit_cooks=True,
        n_cpus=8,
    )

    dds.deseq2()
    stat_res = DeseqStats(dds, n_cpus=8)
    stat_res.summary()
    return stat_res

# Streamlit App Layout
st.title("Differential Gene Expression Analysis")
st.write("""
Welcome to the Differential Gene Expression Analysis page, focusing on Non-Small Cell Lung Cancer (NSCLC) cell lines. This tool allows you to conduct a comprehensive analysis of gene expression alterations in different genetic backgrounds.

## How to Use:
1. **Select a Gene**: Choose a gene from the dropdown menu to perform the DGE analysis on. The gene selection is based on mutations found in the NSCLC cell lines.
2. **Perform Analysis**: Click the 'Perform DGE Analysis' button to start the analysis. This process might take a few moments.
3. **View Results**: After the analysis, the results will be displayed, showing significant gene expression changes. You can also download the full results as a CSV file.

Bellow is an example of the analysis output showing DGE for TP53-mut vs TP53-wt NSCLC cell lines.
""")
# Gene selection dropdown
selected_gene = st.selectbox("Select a Gene for DGE Analysis", mutation_selection)

analysis_performed = False

# Perform DGE analysis button
if st.button('Perform DGE Analysis'):
    analysis_performed = True
    with st.spinner('Running DGE analysis...'):
        dge_results = perform_dge_analysis()
        dge_results = dge_results.results_df
        dge_results = dge_results[dge_results.padj < 0.05]

    # Display the results
    st.write("## DGE Analysis Results")
    st.dataframe(dge_results)

    st.markdown(get_table_download_link(dge_results, "DGE_results.csv", "Download Full DGE Results as CSV"), unsafe_allow_html=True)

# Display the demo file only if the analysis has not been performed
if not analysis_performed:
    st.write("## Demo DGE Results")
    st.dataframe(demo_data)

    st.markdown(get_table_download_link(demo_data, "data/DGE_results_demo.csv", "Download Full DGE Results as CSV"), unsafe_allow_html=True)
