# SAMF
SAMF: Similarity-Aware Multi-omics Fusion

This repository implements the SAMF algorithm — a pipeline for multi-omics data integration and cancer subtype identification.


Module 1: Data Processing

Processing(code1): remove missing values, truncate sample names, remove duplicates

Output:
gene_expression_removal.csv (958 × 59,427)
mirna_removal.csv (962 × 2,234)
TMB_removal.csv (941 × 1)
wei.csv (193 × 16)

Processing(code2): Extracts the common samples across different omics datasets, computes the similarity matrix, and outputs a file as preparation for the self-attention model.


Module 2: Feature Representation via Self-Attention Autoencoder

Training (attention-SNF1.py):
Train a self-attention based autoencoder to learn latent feature representations → saved in models/

Feature Extraction (process_all_data.py):
Use the trained model to extract key feature representations for each omics dataset:
self_attention_omics1.csv
self_attention_omics2.csv
self_attention_omics3.csv
self_attention_omics4.csv

Data Alignment (align.py):
Adjust the sample order of the extracted feature files for downstream analysis:
tmb_aligned.csv
wei_aligned.csv
rna_aligned.csv
dna_aligned.csv

Output: CSV files of key omics feature representations after self-attention processing


Module 3: Clustering & Survival Analysis

Clustering(code3): construct similarity matrix → ICIcluster2.txt

Survival Analysis(code3): input ICIcluster2.txt + clinical.csv → survival.pdf


Module 4: Immune Cell Estimation

Estimation(code4): use CIBERSORT to infer immune cell proportions → CIBERSORT-Results_seed.txt

Visualization: barplot, heatmap, correlation heatmap, violin plot, survival plot


Module 5: Transcription Factor Analysis & logics

Differential Expression Analysis(code5): diffmRNAExp.txt, up.xls, down.xls

Enrichment & Visualization: transcription factor enrichment analysis + heatmaps

Modeling & Evaluation(code6): Apply LASSO regression to select key miRNAs associated with TMB.Build a predictive model using selected features and evaluate performance metrics (AUC).

Outputs: Differentially expressed miRNAs, LASSO feature list, prediction results, and performance statistics.
