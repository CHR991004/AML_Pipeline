# Overview
This repository, AML_pipeline, is dedicated to the comprehensive analysis and modeling of Acute Myeloid Leukemia (AML) using a variety of computational biology techniques. The structure of the repository is designed to systematically address different aspects of AML research, from signature identification to machine learning-based prognostication.

# Directory Structure
The repository is organized into several directories, each focusing on a specific aspect of AML research:

**1. AML signature:** This directory contains resources for the selection of AML signatures (AS), which are crucial for understanding the molecular underpinnings of AML.

**2. RS_DGSE score:** Utilizes rank scores and the DESE algorithm for gene set scoring. This method is instrumental in evaluating the significance of gene expressions in AML.

**3. Enrichment_analysis:** Conducts Gene Ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analyses on AS. This process is vital for understanding the biological functions and pathways associated with AML signatures.

**4. Protein-Protein Interaction (PPI):** Focuses on analyzing the protein-protein interaction networks of AS, providing insights into the complex molecular interactions in AML.

**5. DrugTargetEnrich:** Analyzes potential drug targets within AS through enrichment analysis, aiding in the identification of new therapeutic avenues.

**6. NMF_classification (supplementary):** Applies Non-negative Matrix Factorization (NMF) based classification to AML samples using AS, facilitating the stratification of AML types.

**7. Machine Learning:** Develops a Leave-One-Out Cross-Validation (LOOCV) framework using machine learning algorithms to construct optimal prognostic models for AML patients.

**8. scRNA-seq:** Leverages AS and machine learning algorithms for the classification of AML at the single-cell level, providing a nuanced understanding of AML heterogeneity.

# Purpose
The primary aim of this repository is to serve as a comprehensive resource for researchers and clinicians working on AML. By providing a structured approach to the analysis of AML, it facilitates a deeper understanding of the disease and aids in the development of targeted therapies and prognostic models.

# Usage
Each directory contains input folders and relevant scripts or data files. Users are encouraged to navigate through each directory to understand the specific methodologies and datasets used. The repository is designed to be intuitive for those familiar with computational biology and bioinformatics techniques.

# Contribution and Feedback
Contributions to enhance or extend the capabilities of this repository are welcome. Researchers and developers are encouraged to fork the repository, make improvements, and submit pull requests. For any feedback or queries, please open an issue in the repository.
