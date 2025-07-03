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

### Retrieving sequences using NCBI BLAST.

```bash
python src/blast.py <input_sequence_string> \
                   -o <output_fasta_file> \
                   [-num_aligns <number_of_alignments_to_retrieve>]
```
Example usage:
```bash
python src/blast.py HSQKRVVVLGSGVIGLSSALILARKGYSVHILARDLPEDVSSQTFASPWAGANWTPFMTLTDGPRQAKWEESTFKKWVELVPTGHAMWLKGTRRFAQNEDGLLGHWYKDITPNYRPLPSSECPPGAIGVTYDTLSVHAPKYCQYLARELQKLGATFERRTVTSLEQAFDGADLVVNATGLGAKSIAGIDDQAAEPIRGQTVLVKSPCKRCTMDSSDPASPAYIIPRPGGEVICGGTYGVGDWDLSVNPETVQRILKHCLRLDPTISSDGTIEGIEVLRHNVGLRPARRGGPRVEAERIVLPLDRTKSPLSLGRGSARAAKEKEVTLVHAYGFSSAGYQQSWGAAEDVAQLVDEAFQRYHGAARE -o ./test/results/oxda_rhoto_top_50.fasta -num_aligns 50

```

### Running the Mmseq2 and HMMER workflow.

The main preprocessing workflow is executed via the `run_fido.py` script.

```bash
python run_fido.py -in <input_fasta_file> \
                   -ref <reference_fasta_file> \
                   -out <output_directory> \
                   [-min_seq_id <minimum_sequence_identity_for_clustering>] \
                   [-proj_name <project_name>] \
                   [-add_ref <reference_fasta_to_add_to_reps>]
```
Example usage:
```bash
python run_fido.py -in <input_fasta_file> \
                   -ref <reference_fasta_file> \
                   -out <output_directory> \
                   [-min_seq_id <minimum_sequence_identity_for_clustering>] \
                   [-proj_name <project_name>] \
                   [-add_ref <reference_fasta_to_add_to_reps>]
```

## Configuration

Configuration

The env.py file contains configurable parameters for the preprocessing steps:

    * HMMER_BIN: Path to the HMMER bin directory (e.g., ./hmmer/bin).

    * CLUSTAL_O: Path to the Clustal Omega executable (e.g., ./bin/clustalo).

    * MIN_SEQ_LEN: Minimum sequence length for filtering (default: 0).

    * MAX_SEQ_LEN: Maximum sequence length for filtering (default: 600).

    * MIN_IDEN: Minimum percentage identity for BLAST filtering (default: 30).

    * MAX_IDEN: Maximum percentage identity for BLAST filtering (default: 101, set greater than 100 to include identical sequences).

    * QUERY_COVERAGE: Minimum query coverage for BLAST filtering (default: 50).

These values can be modified in env.py to suit your specific needs. 

## Pipeline Steps

The run_fido.py script executes the following preprocessing steps sequentially:

    1. filter_seqs: Filters input sequences based on length (MIN_SEQ_LEN, MAX_SEQ_LEN) and removes sequences containing 'X' characters. Outputs step_1_filtered.fasta.

    2. blast_filter: Performs a BLASTP search against the filtered sequences using a reference FASTA. It filters sequences based on MIN_IDEN, MAX_IDEN, and QUERY_COVERAGE. Outputs step_2_blastp_filtered.fasta.

    3. run_mmseq: Clusters the BLAST-filtered sequences using MMseqs2 with a specified min_seq_id. Outputs representative sequences to step_3_temp_mmseq_rep_seq.fasta and cluster assignments to step_3_temp_mmseq_cluster.tsv.

    4. Reference Sequence Addition (Optional): If -add_ref is provided, the reference sequence(s) from the specified file are appended to the representative sequences generated by MMseqs2.

    5. clustalo: Performs multiple sequence alignment on the representative sequences (including any added reference sequences) using Clustal Omega. Outputs step_4_rep_seq_aligned.fasta.

    6. hmmer_build_and_align: Builds a Hidden Markov Model (HMM) from the Clustal Omega aligned representative sequences and then aligns the full set of sequences (from step_2_blastp_filtered.fasta) to this HMM. Outputs full_db_alignment.fasta.

    7. build_dataset: Pads the aligned sequences to a uniform length and integrates cluster IDs from the MMseqs2 output. Saves the final dataset as final_alignment_with_clust_ids.csv.

## License

This project is open-source and available under the MIT Open License.

## Contact

For any questions, issues, or contributions, please open an issue on the GitHub repository or contact the maintainers at cookie2004@gmail.com.