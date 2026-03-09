# Anubis

<p align="center">
  <img src="Figure%201.png" width="1000">
</p>

### Ancient Egyptian Collagen Database

**Anubis** is a bioinformatics platform designed for the exploration and analysis of **collagen peptides across species**.  
The project integrates peptide datasets, computational analysis, and an interactive dashboard to study collagen peptide distribution and uniqueness in different organisms.

The platform combines **Python-based analysis pipelines** and an **R Shiny interactive dashboard** to provide researchers with tools for exploring collagen proteomics datasets derived from **UniProt**.

---

# Project Objectives

The main goals of this project are:

- Analyze **collagen peptides across species**
- Compare collagen peptide patterns between **bone and skin**
- Identify **unique peptides per species**
- Build a searchable **collagen peptide database**
- Provide an **interactive visualization dashboard**

This project aims to support bioinformatics research in **proteomics, peptide analysis, and species comparison**.

---

# Data Source

The peptide and protein datasets used in this project are derived from:

**UniProt Protein Database**

The datasets include processed information such as:

- peptide to species mapping
- peptide scoring results
- unique peptide identification
- species summary statistics

These datasets are stored in the `data` directory.

---

# Platform Features

## Peptide Analysis

The Python notebooks perform several computational analyses including:

- peptide scoring
- species mapping
- peptide filtering
- identification of unique peptides
- dataset preprocessing

## Visualization

The Shiny dashboard allows interactive exploration of the data through:

- species distribution visualizations
- peptide counts per species
- protein counts per species
- peptide existence analysis
- species summary statistics

## Data Exploration

Researchers can:

- search for peptides
- explore peptide distribution across species
- compare collagen peptide profiles between **bone and skin**
- visualize proteomics patterns across organisms

---

# Dashboard Components

The **Anubis Shiny Dashboard** contains several analytical sections:

- **Home**
- **Overview**
- **Data Snapshot**
- **Proteins per Species**
- **Peptides per Species**
- **Unique Peptides per Species**
- **Peptide Existence**
- **Existence per Species**
- **Species Summary**
- **Diagnostics**

These modules allow users to explore collagen peptide datasets interactively.

---

# Project Structure



### Directory Description

**data**

Contains processed peptide datasets used for analysis and visualization.

**notebooks**

Jupyter notebooks used for peptide analysis, including:

- peptide scoring
- species mapping
- peptide filtering
- dataset preparation

**shiny**

Contains the **R Shiny dashboard** used for interactive data exploration.

---

# Running the Analysis

## Run the Python Notebooks

The analysis notebooks can be opened in **Jupyter Notebook** or **Google Colab**.

Notebooks available:

- `notebooks/Bone.ipynb`
- `notebooks/Skin.ipynb`

These notebooks perform:

- peptide scoring
- species mapping
- peptide filtering
- dataset preparation

---

## Run the Shiny Dashboard

To launch the interactive dashboard locally:

```r
shiny::runApp("shiny")
Technologies Used

The project combines multiple tools commonly used in bioinformatics:

Python

Jupyter Notebooks

R

R Shiny

Bioinformatics data processing

Proteomics datasets

UniProt protein database

Authors

Hagar Elazab
Department of Biotechnology
Faculty of Agriculture
Cairo University

Nawal Hassan
Undergraduate Studies at the Department of Biotechnology
Faculty of Agriculture
Ain Shams University
Cairo, Egypt
