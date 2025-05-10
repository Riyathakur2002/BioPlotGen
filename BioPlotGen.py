import streamlit as st
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter
import io
import numpy as np

st.set_page_config(page_title="BioPlotGen", layout="wide")
st.title("🧬 BioPlotGen – Bioinformatics Plot Generator")

# ---------------- Navigation Tabs ---------------- #
tab1, tab2, tab3 = st.tabs(["🏠 Welcome", "🧪 About the App", "👥 About Team"])

with tab1:
    st.title("🧬 Welcome to BioPlotGen")
    st.markdown("""
    **BioPlotGen** is an interactive web app for visualizing biological data—whether you're working with gene expression, DNA, or protein sequences. It transforms complex datasets into meaningful plots with just a few clicks.
    """)

with tab2:
    st.header("🔍 About BioPlotGen")
    st.markdown("""
    **BioPlotGen** helps bioinformatics students and researchers generate quick visual insights from biological data.

    ### ⚙️ Features:
    - **📁 File Upload:** Supports .csv, .txt, .fasta
    - **🧬 Sequence Analysis:** GC Content, Codon Usage, Amino Acid Frequency
    - **📊 Plot Types:** Scatter, Line, Bar, Histogram, Box, Violin, Grouped Bar, Heatmap, Volcano
    - **💾 Export Plots:** PNG, PDF, SVG
    """)

with tab3:
    st.header("👩‍💻 About the Creator")
    st.markdown("""
    **Riyasingh Thakur**, MSc Bioinformatics student at DES Pune University.

    *Gratitude to mentors Dr. Kushagra Kashyab and Dr. Poonam Deshpande for their guidance.*

    🔗 *Built using Python, Streamlit, Seaborn, Biopython*
    """)

# ---------------- Helper Functions ---------------- #
def load_csv(file):
    return pd.read_csv(file)

def load_txt(file):
    return pd.read_csv(file, delimiter="\t")

def load_fasta(uploaded_file):
    content = uploaded_file.read().decode("utf-8")
    handle = io.StringIO(content)
    records = list(SeqIO.parse(handle, "fasta"))
    data = [{"ID": rec.id, "Description": rec.description, "Sequence": str(rec.seq)} for rec in records]
    return pd.DataFrame(data)

def calculate_gc(seq):
    seq = seq.upper()
    g = seq.count("G")
    c = seq.count("C")
    return round((g + c) / len(seq) * 100, 2) if len(seq) > 0 else 0

def codon_usage(sequence):
    seq = sequence.upper().replace("\n", "").replace(" ", "")
    codons = [seq[i:i+3] for i in range(0, len(seq) - 2, 3) if len(seq[i:i+3]) == 3]
    return dict(Counter(codons))

def amino_acid_freq(sequence):
    seq_obj = Seq(sequence)
    protein = seq_obj.translate(to_stop=True)
    return dict(Counter(str(protein)))

def volcano_plot(df, x='log2FoldChange', y='pvalue', threshold=0.05, lfc_threshold=1):
    df['Expression'] = np.select(
        [(df[y] < threshold) & (df[x] > lfc_threshold),
         (df[y] < threshold) & (df[x] < -lfc_threshold)],
        ['Overexpressed', 'Underexpressed'],
        default='Not Significant'
    )
    fig, ax = plt.subplots()
    df[y] = df[y].clip(lower=1e-300)  # Prevent log(0)
    sns.scatterplot(data=df, x=x, y=-np.log10(df[y]), hue='Expression', ax=ax,
                    palette={'Overexpressed': 'red', 'Underexpressed': 'blue', 'Not Significant': 'gray'})
    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title("Volcano Plot")
    return fig

def grouped_bar_plot(df, x, y, hue):
    fig, ax = plt.subplots()
    sns.barplot(data=df, x=x, y=y, hue=hue, ax=ax)
    return fig

# ---------------- Sidebar ---------------- #
st.sidebar.header("Upload Your Data")
file = st.sidebar.file_uploader("Upload CSV, TXT, or FASTA", type=["csv", "txt", "fasta"])

plot_types = ["Scatter Plot", "Line Plot", "Bar Chart", "Histogram", "Box Plot", "Violin Plot", "Heatmap", "Grouped Bar Plot", "Volcano Plot"]

# ---------------- File Handling ---------------- #
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

    # ------------- FASTA Handling ------------- #
    if ext == "fasta":
        st.sidebar.header("Sequence Analysis")
        selected_seq_index = st.sidebar.selectbox("Select Sequence", data.index, format_func=lambda x: data.loc[x, 'ID'])
        selected_seq = data.loc[selected_seq_index, 'Sequence']
        analysis_type = st.sidebar.selectbox("Analysis Type", ["Codon Usage", "Amino Acid Frequency"])
        st.subheader(f"{analysis_type} for: {data.loc[selected_seq_index, 'ID']}")

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

    # ------------- Plots for CSV/TXT ------------- #
    else:
        st.sidebar.header("Plot Settings")
        plot_type = st.sidebar.selectbox("Choose Plot Type", plot_types)
        x_col = st.sidebar.selectbox("X-axis", data.columns)
        y_col = st.sidebar.selectbox("Y-axis", data.columns)

        hue_col = None
        if plot_type == "Grouped Bar Plot":
            hue_col = st.sidebar.selectbox("Hue (grouping)", data.columns)

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
                    sns.heatmap(data.select_dtypes(include=['number']).corr(), annot=True, cmap="coolwarm", ax=ax)
                else:
                    st.warning("Heatmap requires at least 2 numeric columns.")
            fig.tight_layout()

        st.pyplot(fig)

        # ---------- Export Plot ---------- #
        st.sidebar.header("Export Plot")
        file_format = st.sidebar.selectbox("Select format", ["PNG", "SVG", "PDF"])
        buf = io.BytesIO()
        fig.savefig(buf, format=file_format.lower(), bbox_inches="tight")
        st.download_button(
            label="Download Plot",
            data=buf.getvalue(),
            file_name=f"bioplotgen_plot.{file_format.lower()}",
            mime=f"image/{file_format.lower()}"
        )

