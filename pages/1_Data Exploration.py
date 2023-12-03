import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

st.title('Data Exploration')

# Introductory text about the page
st.write("""
Explore various genetic data aspects through interactive visualizations. Customize your view by making selections in the sidebar and dive into the analysis of Oncotree Primary Diseases and common gene mutations.
""")


# Load necessary data
model = pd.read_csv('data/Model.csv')
mutation = pd.read_csv('data/OmicsSomaticMutations_filter.csv', index_col=0)

# First Plot: Top N OncotreePrimaryDisease Counts
st.subheader("OncotreePrimaryDisease Counts")
st.write("""
View the counts of different Oncotree Primary Diseases. Adjust the number of top diseases to display using the slider in the sidebar. This feature is useful for identifying the most prevalent diseases in the dataset.
""")


top_n = st.sidebar.slider("Select Top N to Show", min_value=1, max_value=20, value=10)
OncotreePrimaryDisease_count = model.OncotreePrimaryDisease.value_counts()
plt.figure()
plt.barh(OncotreePrimaryDisease_count.head(top_n).index, OncotreePrimaryDisease_count.head(top_n))
for index, value in enumerate(OncotreePrimaryDisease_count.head(top_n)):
    plt.text(value, index, str(value))
plt.xlabel('Count')
plt.ylabel('')
plt.title(f'Top {top_n} OncotreePrimaryDisease Counts')
st.pyplot(plt)

# Second Plot: Most Common Mutated Genes
st.subheader("Most Common Mutated Genes")
st.write("""
Explore the most frequently mutated genes. Choose a specific Oncotree Primary Disease from the sidebar to narrow down the analysis, or select 'All' to view the overall mutation landscape. This visualization helps in understanding gene mutation patterns across different diseases.
""")

disease_options = list(model['OncotreePrimaryDisease'].unique())
disease_options.insert(0, 'All')
disease_selection = st.sidebar.selectbox('Select OncotreePrimaryDisease', disease_options)

if disease_selection != 'All':
    model_filter = model[model['OncotreePrimaryDisease'] == disease_selection]
    mutation_filter = mutation[mutation['ModelID'].isin(model_filter['ModelID'])]
else:
    mutation_filter = mutation

mutation_filter = mutation_filter[(mutation_filter.Driver == True) | (mutation_filter.LikelyDriver == True)]

plt.figure()
plt.barh(mutation_filter.HugoSymbol.value_counts().head(top_n).index, mutation_filter.HugoSymbol.value_counts().head(top_n))
for index, value in enumerate(mutation_filter.HugoSymbol.value_counts().head(top_n)):
    plt.text(value, index, str(value))
plt.xlabel('Count')
plt.ylabel('')
title = 'Top {} Most Common Mutated Genes'.format(top_n)
if disease_selection != 'All':
    title += f' ({disease_selection})'
plt.title(title)
st.pyplot(plt)


  
