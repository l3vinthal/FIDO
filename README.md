# FIDO: Protein Sequence Preprocessing for Deep Learning

FIDO is a suite of scripts designed to preprocess protein sequences for downstream deep learning applications. It automates common bioinformatics tasks such as sequence filtering, BLAST-based filtering, sequence clustering using MMseqs2, multiple sequence alignment using Clustal Omega, and HMMER-based alignment.

## Table of Contents

- [FIDO: Protein Sequence Preprocessing for Deep Learning](#fido-protein-sequence-preprocessing-for-deep-learning)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Installation](#installation)
    - [Prerequisites](#prerequisites)
    - [Setup Steps](#setup-steps)
  - [Usage](#usage)
    - [Running the Preprocessing Workflow](#running-the-preprocessing-workflow)
    - [Example Notebook](#example-notebook)
  - [Configuration](#configuration)
  - [Pipeline Steps](#pipeline-steps)
  - [License](#license)
  - [Contact](#contact)

## Features

* **Sequence Filtering**: Filters sequences based on length and presence of 'X' characters.
* **BLAST-based Filtering**: Filters sequences based on identity and query coverage against a reference sequence.
* **MMseqs2 Clustering**: Clusters sequences to reduce redundancy.
* **Clustal Omega Alignment**: Performs multiple sequence alignment on representative sequences from clustering.
* **HMMER Alignment**: Builds an HMM from representative sequences and aligns the full database to it.
* **Dataset Generation**: Creates a CSV file with aligned sequences and their cluster IDs, padded for consistent length.

## Installation

### Prerequisites

* `conda` (Miniconda or Anaconda) is required to manage the environment and dependencies.

### Setup Steps

1.  **Clone the Repository**:
    ```bash
    git clone github.com/l3vinthal/FIDO
    cd FIDO
    ```

2.  **Run the Environment Setup Script**:
    The `env_setup.sh` script automates the creation of a `conda` environment named `fido` and installs necessary bioinformatics tools like HMMER and Clustal Omega, along with Python packages.

    ```bash
    bash env_setup.sh
    ```

    This script will perform the following actions:
    * Create an installation directory `$HOME/`.
    * Create a `conda` environment named `fido` with Python 3.10.
    * Install HMMER (version 3.4)
    * Install Clustal Omega (version 1.2.4)
    * Install `blast` and `mmseqs2` via Bioconda.
    * Install Python packages like `pandas` via conda.
    * Clean up downloaded archives.

3.  **Activate the Conda Environment**:
    After installation, activate the `fido` environment:
    ```bash
    conda activate fido
    ```

## Usage

### Running the Preprocessing Workflow

The main preprocessing workflow is executed via the `run_fido.py` script.

```bash
python run_fido.py -in <input_fasta_file> \
                   -ref <reference_fasta_file> \
                   -out <output_directory> \
                   [-min_seq_id <minimum_sequence_identity_for_clustering>] \
                   [-proj_name <project_name>] \
                   [-add_ref <reference_fasta_to_add_to_reps>]