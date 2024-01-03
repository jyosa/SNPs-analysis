# SNPs-analysis


## Objective
This script is designed to analyze amino acid mutations from MSA data. It focuses on categorizing amino acids based on their chemical properties, performing clustering, and generating visualizations and CSV files for further analysis.

## Source of MSA Data
The MSA data used in this analysis is obtained from the Anthony Nolan GitHub repository, which can be accessed [here](https://github.com/ANHIG/IMGTHLA). The data is in a specific format different from standard alignment formats like Clustal or Muscle, making it incompatible with libraries like Biopython. This tool is developed to handle this unique format.

## Algorithm Description

### File Reading and Initial Processing
- Reads an MSA file.
- Extracts protein blocks and converts them into a format suitable for analysis.

### Data Preparation and Transformation
- Concatenates sequences from different allele blocks.
- Fills sequences with zeros to match the maximum length.
- Numbers the amino acids in the reference sequence.
- Finds the index of a specific amino acid and slices the data accordingly.

### Data Analysis
- Transforms the data based on certain rules (handling of '-', '.', '*', ‘X'). Rules are in line with [EBI standards](https://www.ebi.ac.uk/ipd/imgt/hla/alignment/help/).
- Categorizes amino acids by their chemical properties.
- Creates a DataFrame for analysis, including frequency distributions and clustering using k-modes.

### Visualization
- Generates histograms and interactive plots to visualize the data analysis results.

### Cluster Analysis and File Output
- Saves each cluster's data into separate CSV files for further analysis.

### Command Line Interface (CLI)
- Provides a CLI for users to specify input parameters.

## Inputs
- **MSA file**: A file containing multiple sequence alignments.
- **aa_init**: Initial amino acid position in the reference sequence.
- **aa_pos**: Position of amino acids for analysis.

## Outputs

### Visualizations
- Histogram of relative frequency percentages of amino acids.
- Stacked bar plot showing the distribution of chemical properties in each cluster.

### CSV Files
- A file named 'aminoacid_list.csv' containing categorized amino acids and their chemical properties, lists of amino acids, and their alleles.
- Separate CSV files for each cluster with detailed data.

### Terminal Output
- Printout of the DataFrame with cluster labels.
- Instructions and usage information (if requested).

## Detailed Functions
- `extract_prot_blocks_with_format`: Extracts protein blocks from file content.
- `convert_blocks_to_dict_filtered`: Converts protein blocks to a dictionary, filtering specific parts.
- `concatenate_sequences_from_dict_list`: Concatenates allele sequences from a list of dictionaries.
- `add_position_numbers_to_reference`: Adds position numbers to the reference sequence.
- `get_indx`: Finds the index of an amino acid in a sequence.
- `fill_with_zeros`: Fills a string with zeros to reach a desired length.
- `number_letters`: Numbers letters in a list starting from a specified number.
- `find_index_by_number`: Finds the index of an item in a list based on its number.
- `slice_dictionary_at_index`: Slices a dictionary at a specific index.
- `transform_value`: Transforms a value based on specific rules.
- `categorize_amino_acid`: Categorizes an amino acid based on its type.
- `process`: The main function for processing and analysis.
- `print_help`: Prints help and usage instructions.
- `main`: The entry point of the script, handling CLI.

## Downloading the Input MSA Data
The input MSA data can be downloaded from the [following link](https://github.com/ANHIG/IMGTHLA/blob/Latest/Alignments_Rel_3540.zip).

# Installation and Usage Guide

## Prerequisites
Before running the script, ensure you have Python installed on your system. Python 3.6 or newer is recommended. You can download and install Python from [python.org](https://www.python.org/downloads/).

## Required Libraries
The script requires several Python libraries. These can be installed using `pip`, the Python package manager. The required libraries are:
- pandas
- matplotlib
- scipy
- scikit-learn
- kmodes
- plotly

You can install these libraries by running the following command:

```bash
pip install pandas matplotlib scipy scikit-learn kmodes plotly
```

## Downloading the Script

Download the script from the GitHub repository or clone it using the following command:

```bash
git clone https://github.com/jyosa/SNPs-analysis.git
```

## Running the Script

To run the script, navigate to the directory containing the script and run it using Python. The script requires three command-line arguments:

    -f: Path to the MSA file.
    -ai: Initial amino acid position in the reference sequence.
    -ap: Position of the amino acids for analysis.

For example:

The reference sequence from HLA_B starts at -23, which includes the signal peptide, so we keep that to follow the amino acid indexation 

```bash
python hla2.py -f 'B_prot.txt' -ai -23 -ap 116
```

## Mutation at 116 position

![Fig 1. SNPs distribution and cluster by chemical properties for position 116]([relative/path/in/repository/to/image.svg](https://github.com/jyosa/SNPs-analysis/blob/main/img.001.png))

Regarding the significance of position 116 in HLA-B according to results: This gene encodes for a key component of the Class I major histocompatibility complex (MHC), essential in immune responses. HLA-B molecules have specific peptide-binding motifs determined by the structure of their peptide-binding pockets, known as the B and F pockets. Position 116 and other amino acids within these pockets are crucial for determining the binding specificity of peptides. These peptides are then presented to T cells, a process central to the immune system's ability to recognize and respond to pathogens. The specific amino acids at positions like 116 influence the binding and presentation of these peptides, thereby playing a pivotal role in modulating immune responses.

Mutations at position 116, although less common, are of particular interest. For example, mutations leading to a stop codon or the incorporation of proline can significantly alter the cleft-like structure of the peptide-binding region. This alteration can impact the overall function of the HLA-B molecule. Similarly, the substitution of arginine at this position can affect the affinity of the peptide-binding pocket for its ligands. These changes are crucial in understanding variations in immune response.

In the context of autoimmune diseases, mutations at position 116 and their correlation with specific peptides presented to T cells can be particularly enlightening. Exploring these mutations further is essential, especially in association studies that correlate such genetic variations (SNPs) with autoimmune disorders. Such research could provide valuable insights into the genetic underpinnings of these diseases and open avenues for targeted therapies.


Powered by @El_Dryosa
