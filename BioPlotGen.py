import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import io
import numpy as np

# Setting page configuration
st.set_page_config(page_title="BioPlotGen", layout="wide")

# Apply light pastel background and custom styles
st.markdown("""
    <style>
    body {
        background-color: #f4f9f9;
    }
    .stApp {
        background-color: #fefefe;
        color: #333333;
    }
    .css-18e3th9 {
        padding: 2rem;
        border-radius: 1rem;
        background: #f0f8ff;
    }
    h1, h2, h3 {
        color: #005f73;
        font-weight: bold;
    }
    h1 {
        font-size: 40px;
        color: #ff6347;  # Tomato color for title
    }
    .stButton>button {
        background-color: #94d2bd;
        color: black;
        border-radius: 10px;
        padding: 10px 16px;
    }
    .stDownloadButton>button {
        background-color: #ffb703;
        color: black;
        border-radius: 10px;
        padding: 10px 16px;
    }
    .stSelectbox>div>div {
        background-color: #ffffff;
    }
    .stTabs>div {
        background-color: #e0f7fa;
        color: #004d40;
    }
    </style>
""", unsafe_allow_html=True)

# ---------------- Navigation Tabs ---------------- #

# Define tabs
tabs = st.selectbox("Select Tab", ["Welcome", "About the App", "About the Creator", "Plot Generator"])

# Tab 1: Welcome
if tabs == "Welcome":
    st.title("🧬 Welcome to BioPlotGen")
    st.markdown("""
    Welcome to BioPlotGen, the ultimate tool for visualizing and analyzing biological data, designed for biologists, researchers, and students at all levels—no coding required!
    BioPlotGen makes it easier than ever to transform complex biological data into insightful, publication-quality plots with just a few clicks. Whether you're analyzing gene expression, DNA sequences, protein structures, or other bioinformatics data, BioPlotGen provides an intuitive, user-friendly platform to help you make sense of your research without needing any programming skills.
    """)

# Tab 2: About the App
elif tabs == "About the App":
    st.header("🔍 About BioPlotGen")
    st.markdown("""
    **BioPlotGen** helps bioinformatics students and researchers generate quick visual insights from biological data.

    ### ⚙️ Features:
    - **📁 File Upload:** Supports .csv, .txt, .fasta
    - **🧬 Sequence Analysis:** GC Content, Codon Usage, Amino Acid Frequency
    - **📊 Plot Types:** Scatter, Line, Bar, Histogram, Box, Violin, Grouped Bar, Heatmap, Volcano
    - **💾 Export Plots:** PNG, PDF, SVG
    """)

# Tab 3: About the Creator
elif tabs == "About the Creator":
    st.header("👩‍💻 About the Creator")
    st.markdown("""
    **Riyasingh Thakur**, MSc Bioinformatics student at DES Pune University.
    I am a passionate MSc Bioinformatics student at DES Fergusson College, Pune, with a strong background in microbiology and a growing expertise in programming. I am deeply committed to bridging the gap between biological sciences and computational technologies.
    My work on BioPlotGen reflects my enthusiasm for creating user-friendly tools that simplify complex biological data analysis. I believe that the integration of biology and computer science holds the key to unlocking breakthroughs in modern research, and I am determined to be a part of that transformation.

    I would like to express my deepest gratitude to Dr. Kushagra Kashyab and Dr. Poonam Deshpande for their invaluable guidance, encouragement, and continuous support throughout the development of BioPlotGen.
    Their insights in the field of bioinformatics and constant motivation have been instrumental in shaping this project from concept to execution.
    I also extend sincere thanks to my peers, faculty, and the Department of Bioinformatics, DES Fergusson College, Pune, for providing a collaborative environment and helpful feedback at every stage of development.
    Finally, I am thankful to the open-source community and tools like Python, Streamlit, Seaborn, Biopython, and Matplotlib, which made the creation of this interactive application possible
    🔗 *Built using Python, Streamlit, Seaborn, Biopython*
    """)

# Tab 4: Plot Generator
elif tabs == "Plot Generator":
    st.header("📈 Plot Generator")

    file = st.file_uploader("Upload CSV, TXT, or FASTA", type=["csv", "txt", "fasta"])

    # Function Definitions

    def load_csv(file):
        return pd.read_csv(file)

    def load_txt(file):
        return pd.read_csv(file, sep="\t")  # or sep="\s+" for space-separated

    def load_fasta(file):
        records = []
        for record in SeqIO.parse(file, "fasta"):
            records.append({"ID": record.id, "Sequence": str(record.seq)})
        return pd.DataFrame(records)

    def calculate_gc(sequence):
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100 if sequence else 0

    def codon_usage(sequence):
        codons = [sequence[i:i+3] for i in range(0, len(sequence)-2, 3)]
        return dict(Counter(codons))

    def amino_acid_freq(sequence):
        protein_seq = str(Seq(sequence).translate(to_stop=True))
        return dict(Counter(protein_seq))

    def grouped_bar_plot(df, x, y, hue):
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.barplot(data=df, x=x, y=y, hue=hue, ax=ax)
        return fig

    def volcano_plot(df, logFC_col, pval_col):
        df['-log10(pval)'] = -np.log10(df[pval_col])
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.scatterplot(data=df, x=logFC_col, y='-log10(pval)', hue=(df[pval_col] < 0.05) & (abs(df[logFC_col]) > 1), palette={True: 'red', False: 'gray'}, ax=ax)
        ax.axhline(-np.log10(0.05), color='blue', linestyle='--')
        ax.axvline(-1, color='blue', linestyle='--')
        ax.axvline(1, color='blue', linestyle='--')
        return fig

    plot_types = ["Scatter Plot", "Line Plot", "Bar Chart", "Histogram", "Box Plot", "Violin Plot", "Grouped Bar Plot", "Heatmap", "Volcano Plot"]

    if file:
        ext = file.name.split(".")[-1].lower()
        if ext == "csv":
            data = load_csv(file)
        elif ext == "txt":
            data = load_txt(file)
        elif ext == "fasta":
            data = load_fasta(file)
            data['GC_Content (%)'] = data['Sequence'].apply(calculate_gc)

        st.subheader("📄 Uploaded Data Preview")
        st.dataframe(data.head())

        if ext == "fasta":
            st.subheader("🧬 Sequence Analysis")
            selected_seq_index = st.selectbox("Select Sequence", data.index, format_func=lambda x: data.loc[x, 'ID'])
            selected_seq = data.loc[selected_seq_index, 'Sequence']
            analysis_type = st.selectbox("Analysis Type", ["Codon Usage", "Amino Acid Frequency"])
            st.write(f"**{analysis_type} for:** `{data.loc[selected_seq_index, 'ID']}`")

            if analysis_type == "Codon Usage":
                codon_freq = codon_usage(selected_seq)
                codon_df = pd.DataFrame(codon_freq.items(), columns=["Codon", "Count"]).sort_values(by="Count", ascending=False)
                st.dataframe(codon_df)
                fig, ax = plt.subplots(figsize=(12, 6))
                sns.barplot(data=codon_df, x="Codon", y="Count", palette="viridis", ax=ax)
                plt.xticks(rotation=45)
                st.pyplot(fig)

            elif analysis_type == "Amino Acid Frequency":
                aa_freq = amino_acid_freq(selected_seq)
                aa_df = pd.DataFrame(aa_freq.items(), columns=["Amino Acid", "Count"]).sort_values(by="Count", ascending=False)
                st.dataframe(aa_df)
                fig, ax = plt.subplots(figsize=(12, 6))
                sns.barplot(data=aa_df, x="Amino Acid", y="Count", palette="plasma", ax=ax)
                plt.xticks(rotation=45)
                st.pyplot(fig)

        else:
            st.subheader("📊 Plot Settings")
            plot_type = st.selectbox("Choose Plot Type", plot_types)
            x_col = st.selectbox("X-axis", data.columns)
            y_col = st.selectbox("Y-axis", data.columns)

            hue_col = None
            if plot_type == "Grouped Bar Plot":
                hue_col = st.selectbox("Hue (grouping)", data.columns)

            st.subheader(f"{plot_type} of {x_col} vs {y_col}")

            if plot_type == "Grouped Bar Plot":
                fig = grouped_bar_plot(data, x_col, y_col, hue_col)
            elif plot_type == "Volcano Plot":
                fig = volcano_plot(data, x_col, y_col)
            else:
                fig, ax = plt.subplots()
                if plot_type == "Scatter Plot":
                    sns.scatterplot(data=data, x=x_col, y=y_col, ax=ax)
                elif plot_type == "Line Plot":
                    sns.lineplot(data=data, x=x_col, y=y_col, ax=ax)
                elif plot_type == "Bar Chart":
                    sns.barplot(data=data, x=x_col, y=y_col, ax=ax)
                elif plot_type == "Histogram":
                    sns.histplot(data[x_col], kde=True, ax=ax)
                elif plot_type == "Box Plot":
                    sns.boxplot(data=data, x=x_col, y=y_col, ax=ax)
                elif plot_type == "Violin Plot":
                    sns.violinplot(data=data, x=x_col, y=y_col, ax=ax)
                elif plot_type == "Heatmap":
                    if data.select_dtypes(include=['number']).shape[1] >= 2:
                        corr_matrix = data.corr()
                        fig, ax = plt.subplots(figsize=(10, 8))
                        sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", fmt=".2f", ax=ax)
                        st.pyplot(fig)
                st.pyplot(fig)

            st.subheader("💾 Download Plot")
            st.download_button("Download Plot", data=fig, file_name="plot.png", mime="image/png")

        st.subheader("🚀 Try another File!")
        st.button("Reset")


          
