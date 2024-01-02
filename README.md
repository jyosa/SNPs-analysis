# SNPs-analysis
# Amino Acid Sequence Analysis Tool

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
