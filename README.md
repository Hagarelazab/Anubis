# Anubis
![Python](https://img.shields.io/badge/Python-Analysis-blue)
![R](https://img.shields.io/badge/R-Shiny-green)
![Bioinformatics](https://img.shields.io/badge/Field-Bioinformatics-orange)
![License](https://img.shields.io/badge/License-GPLv3-red)
<p align="center">
  <img src="./Figure%201.png" width="1000">
</p>

# Anubis - Ancient Egyptian Collagen Database

**Anubis** is a bioinformatics platform designed for the exploration and analysis of **collagen peptides across species**.

The project integrates peptide datasets, computational analysis, and an interactive **R Shiny dashboard** to study collagen peptide distribution, peptide uniqueness, and species-level collagen patterns using data derived from the **UniProt Protein Database**.

![Anubis Workflow](Figure%201.png)

**Figure 1.** Workflow of the Anubis platform showing data collection from UniProt, in-silico peptide digestion, peptide scoring, and interactive visualization in the Shiny dashboard.

---

## Project Objectives

The main goals of this project are to:

* Analyze collagen peptides across species.
* Compare collagen peptide patterns between bone and skin.
* Identify unique peptides per species.
* Build a searchable collagen peptide database.
* Provide an interactive visualization dashboard for collagen proteomics exploration.

This project supports bioinformatics research in **proteomics, peptide analysis, collagen biology, and species comparison**.

---

## Data Source

The peptide and protein datasets used in this project are derived from:

**UniProt Protein Database**

The processed datasets include:

* Peptide-to-species mapping.
* Peptide scoring results.
* Unique peptide identification.
* Species summary statistics.
* Bone and skin collagen peptide datasets.

These datasets are stored in the `data/` directory.

---

## Platform Features

### Peptide Analysis

The Python notebooks perform several computational analyses, including:

* Peptide scoring.
* Species mapping.
* Peptide filtering.
* Identification of unique peptides.
* Dataset preprocessing and preparation.

### Interactive Visualization

The Shiny dashboard allows users to explore the processed datasets through:

* Species distribution visualizations.
* Peptide counts per species.
* Protein counts per species.
* Unique peptide summaries.
* Peptide existence analysis.
* Species-level summary statistics.

### Data Exploration

Researchers can use the platform to:

* Search for peptides.
* Explore peptide distribution across species.
* Compare collagen peptide profiles between bone and skin.
* Visualize proteomics patterns across organisms.
* Investigate species-specific collagen peptide uniqueness.

---

## Dashboard Components

The **Anubis Shiny Dashboard** contains several analytical sections:

* Home
* Overview
* Data Snapshot
* Proteins per Species
* Peptides per Species
* Unique Peptides per Species
* Peptide Existence
* Existence per Species
* Species Summary
* Diagnostics

These modules allow users to explore collagen peptide datasets interactively.

---

## Project Structure

```text
ANUBIS/
├── data/
│   ├── step2_peptide_species_map.csv
│   ├── step2_peptides_with_species.csv
│   ├── step3_species_summary.csv
│   └── step3_unique_peptides_scored.csv
│
├── notebooks/
│   ├── Bone.ipynb
│   └── Skin.ipynb
│
├── shiny/
│   └── app.R
│
├── Figure 1.png
├── LICENSE
└── README.md
```

---

## Directory Description

### `data/`

Contains the processed peptide datasets used for analysis and visualization.

### `notebooks/`

Contains the Jupyter notebooks used for peptide analysis, including:

* Peptide scoring.
* Species mapping.
* Peptide filtering.
* Dataset preparation.

### `shiny/`

Contains the **R Shiny dashboard** used for interactive data exploration.

---

## Running the Analysis

### Run the Python Notebooks

The analysis notebooks can be opened in **Jupyter Notebook** or **Google Colab**.

Available notebooks:

* `notebooks/Bone.ipynb`
* `notebooks/Skin.ipynb`

These notebooks perform peptide scoring, species mapping, peptide filtering, and dataset preparation.

---

## Run the Shiny Dashboard

To launch the interactive dashboard locally, first install the required R packages:

```r
install.packages(c(
  "shiny", "bslib", "dplyr", "readr", "DT", "plotly", "ggplot2",
  "stringr", "purrr", "scales", "tidyr", "htmltools", "tibble", "xml2"
))
```

Then run the Shiny app:

```r
shiny::runApp("shiny")
```

The required CSV/XML data files should be placed in the `data/` folder.

---

## Technologies Used

This project combines multiple tools commonly used in bioinformatics and proteomics research:

* Python
* Jupyter Notebook
* R
* R Shiny
* UniProt Protein Database
* Bioinformatics data processing
* Proteomics dataset analysis
* Interactive data visualization

---

## Authors

**Hagar Elazab**
Department of Biotechnology
Faculty of Agriculture
Cairo University

**Nawal Hassan**
Undergraduate Studies at the Department of Biotechnology
Faculty of Agriculture
Ain Shams University
Cairo, Egypt

---

## License

This project is licensed under the terms of the license included in this repository.

