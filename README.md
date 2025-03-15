
README file

This repository contains a Shiny application and a set of scripts to process and analyze population genetics data, specifically focusing on SNP (Single Nucleotide Polymorphism) data for detecting Runs of Homozygosity (ROH). The app performs statistical analysis, visualizations, and allows comparisons of genetic datasets.

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Required Libraries](#required-libraries)
4. [File Processing](#file-processing)
5. [Usage](#usage)
6. [License](#license)

## Introduction
**ShinyROH** is an interactive web application designed to process and analyze SNP data to detect ROH (Runs of Homozygosity). It calculates ROH lengths and LOD scores to classify different ROH segments. By leveraging parallel computing, ShinyROH efficiently processes large SNP datasets.

The app allows users to:

- Upload SNP genotype data in the given format:

| rsid          | chromosome | position | allele1 | allele2 |
|---------------|------------|----------|---------|---------|
| rs547237130   | 1          | 72526    | A       | A       |
| rs562180473   | 1          | 565703   | A       | A       |
| rs575203260   | 1          | 567693   | T       | T       |
| rs3131972     | 1          | 752721   | A       | G       |
| rs200599638   | 1          | 752918   | G       | G       |

- Process data using bootstrap-based allele frequency estimation, which adds an extra column (`genotype`) to the table:

| rsid          | chromosome | position | allele1 | allele2 | genotype |
|---------------|------------|----------|---------|---------|----------|
| rs547237130   | 1          | 72526    | A       | A       | 0        |
| rs562180473   | 1          | 565703   | A       | A       | 0        |
| rs575203260   | 1          | 567693   | T       | T       | 0        |
| rs3131972     | 1          | 752721   | A       | G       | 1        |
| rs200599638   | 1          | 752918   | G       | G       | 0        |


- Compute ROH classifications using a Gaussian Mixture Model (GMM)
- Visualize ROH data using interactive plots
- Compare ROH lengths and distributions across chromosomes
- Place user input data in a worldwide ROH distribution plot

## Installation

### Prerequisites
Before running the app, ensure that you have the following software installed:

- **R** version 4.4.2 or higher
- **RStudio** (latest version recommended)

### Required libraries
To run this app, ensure that you have the following R libraries installed:

```r
install.packages(c("shiny", "data.table", "gtools", "tidyr", "dplyr", "mclust", 
                   "MASS", "ggplot2", "plotly", "DT", "future.apply"))


