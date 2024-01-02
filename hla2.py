from itertools import islice
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import OneHotEncoder
from kmodes.kmodes import KModes
import plotly.express as px
import sys

# Function to extract protein blocks from a file's content
def extract_prot_blocks_with_format(file_content):
    lines = file_content.split('\n')  # Splitting file content into lines
    prot_blocks = []  # List to store protein blocks
    current_block = []  # Temporary list to hold current block

    # Loop through each line in the file
    for line in lines:
        if 'Prot' in line:  # Check if line contains 'Prot'
            if current_block:
                # Add the current block to prot_blocks and reset it
                prot_blocks.append('\n'.join(current_block))
                current_block = []
        # If it's part of a block or the start of a new one, add the line to current_block
        if current_block or 'Prot' in line:
            current_block.append(line)

    # Add the last block if it exists
    if current_block:
        prot_blocks.append('\n'.join(current_block))

    return prot_blocks

# Function to convert the protein blocks to a dictionary, filtered based on certain criteria
def convert_blocks_to_dict_filtered(prot_blocks):
    allele_sequences = {}  # Dictionary to hold allele sequences

    # Process each block
    for block in prot_blocks:
        lines = block.split('\n')
        # Iterate starting from the third line in each block
        for line in lines[2:]:
            parts = line.split(maxsplit=1)
            if len(parts) < 2 or parts[0].startswith('Please'):
                continue

            # Split line into allele and sequence, then add to the dictionary
            allele, sequence = parts[0], parts[1]
            allele_sequences[allele] = sequence

    return allele_sequences

# Function to concatenate sequences from a list of dictionaries
def concatenate_sequences_from_dict_list(dict_list):
    concatenated_sequences = {}  # Dictionary to hold concatenated sequences

    # Process each dictionary in the list
    for allele_dict in dict_list:
        for allele, sequences in allele_dict.items():
            # Concatenate sequences or initialize them in the dictionary
            if allele not in concatenated_sequences:
                concatenated_sequences[allele] = sequences
            else:
                concatenated_sequences[allele] += " " + sequences

    return concatenated_sequences

# Add position numbers to the reference sequence in the dictionary
def add_position_numbers_to_reference(sequences, start_position):
    reference_allele = next(iter(sequences))  # Get the first allele as the reference
    reference_sequence = sequences[reference_allele]  # Get the reference sequence

    position = start_position  # Initialize position
    numbered_sequence = ""  # String to hold numbered sequence

    # Loop through each character in the sequence
    for char in reference_sequence:
        if char.isalpha():
            # If character is a letter, append it with its position
            numbered_sequence += f"{char}{position} "
            position += 1
        else:
            numbered_sequence += f"{char} "  # Append non-letter characters as is

    # Update the sequence in the dictionary with numbered sequence
    sequences[reference_allele] = numbered_sequence.strip()

    return sequences

# Function to find the index of an amino acid in a sequence
def get_indx(dict_seq, ind_aa):
    first_key, first_value = next(iter(dict_seq.items()))  # Get the first sequence
    seq_val_list = first_value.split(' ')  # Split the sequence into a list

    # Loop through the list to find the amino acid
    for item_list in seq_val_list:
        if len(item_list) > 2 and int(item_list[1:]) == ind_aa: 
            return seq_val_list.index(item_list)

# Function to fill a string with zeros to reach a desired length
def fill_with_zeros(s, desired_length):
    s = s.replace(' ', '0')  # Replace spaces with zeros
    zeros_needed = desired_length - len(s)  # Calculate needed zeros

    # Fill the string with zeros
    if zeros_needed > 0:
        s += '0' * zeros_needed
    return s

# Function to number the letters in a list starting from a specified number
def number_letters(input_list, start_num):
    new_list = []
    counter = start_num  # Initialize counter

    # Loop through the elements of the list
    for element in input_list:
        if element.isalpha():  # Check if element is a letter
            new_list.append(f'{element}{counter}')  # Append letter with counter
            counter += 1
        else:
            new_list.append(element)  # Append non-letter elements as is

    return new_list

# Function to find the index of an item in a list based on its number
def find_index_by_number(lst, number):
    for i, item in enumerate(lst):
        num_part = ''.join(filter(str.isdigit, item))  # Extract numerical part
        if num_part.isdigit() and int(num_part) == number:
            return i
    return "Number not found in the list."

# Function to slice a dictionary at a specific index
def slice_dictionary_at_index(original_dict, index):
    new_dict = {}
    for key, value_list in original_dict.items():
        # Slice the list at the given index, or set None if index is out of range
        new_dict[key] = value_list[index] if len(value_list) > index else None
    return new_dict

# Function to transform a value based on specific rules
def transform_value(val, original):
    if val == '-':
        return original
    elif val == '.':
        return 'INDEL'
    elif val == '*':
        return 'unknown'
    elif val == 'X':
        return 'Stop'
    else:
        return val


def categorize_amino_acid(aa):
    if aa == 'X' or aa == '0':  # Handle unknown amino acids
        return 'Unknown'

    # Define categories
    aliphatic = ['G', 'A', 'V', 'L', 'I']
    aromatic = ['F', 'Y', 'W']
    acidic = ['D', 'E']
    basic = ['K', 'R', 'H']
    sulfur_containing = ['C', 'M']
    alcohol_containing = ['S', 'T']
    amide_containing = ['N', 'Q']
    imidazole_containing = ['H']
    guandino_containing = ['R']

    # Categorization
    if aa in aliphatic:
        return 'Aliphatic'
    elif aa in aromatic:
        return 'Aromatic'
    elif aa in acidic:
        return 'Acidic'
    elif aa in basic:
        return 'Basic'
    elif aa in sulfur_containing:
        return 'Sulfur-containing'
    elif aa in alcohol_containing:
        return 'Alcohol-containing'
    elif aa in amide_containing:
        return 'Amide-containing'
    elif aa in imidazole_containing:
        return 'Imidazole-containing'
    elif aa in guandino_containing:
        return 'Guandino-containing'
    else:
        return 'Other'

def process(msa,aa_init,aa_pos):

    file_path = msa
    with open(file_path, 'r') as file:
        file_content = file.read()


    # Initial setup and processing
    aa_init = aa_init
    prot_blocks_with_format = extract_prot_blocks_with_format(file_content)
    allele_sequence_dicts_filtered = [convert_blocks_to_dict_filtered([block]) for block in prot_blocks_with_format]
    concatenated_sequences = concatenate_sequences_from_dict_list(allele_sequence_dicts_filtered)

    # Finding the maximum sequence length
    max_num = [len(value) for key, value in concatenated_sequences.items()]
    max_num_t = max(max_num)
    

    # Filling sequences with zeros to match the maximum length
    new_dict = {key2: list(fill_with_zeros(value2, max_num_t)) for key2, value2 in concatenated_sequences.items()}

    # Numbering the letters in the reference sequence
    init_num = aa_init
    reference = next(iter(new_dict.values()))
    refer_list_mark = number_letters(list(reference), init_num)

    # Finding a specific amino acid position
    aa_pos = aa_pos
    index_position = find_index_by_number(refer_list_mark, aa_pos)

    # Slicing the dictionary at the found index
    final_data = slice_dictionary_at_index(new_dict, index_position)
   

    # Converting the final data to a DataFrame 
    df_final_data = pd.DataFrame(list(final_data.items()), columns=['Allele', 'AA'])

    # Applying the transformation to the DataFrame
    original_letter = df_final_data['AA'].iloc[0]
    df_final_data['AA'] = df_final_data['AA'].apply(lambda x: transform_value(x, original_letter))
    

    # Apply the categorization function to the DataFrame
    df_final_data['Chemical Property'] = df_final_data['AA'].apply(categorize_amino_acid)
    df_final_data.to_csv('aminoacid_list.csv', index=False)
    # Plotting a histogram of frequencies for the values in column 2
    value_counts = df_final_data['AA'].value_counts(normalize=True) * 100
    value_counts.plot(kind='bar')
    plt.title('Percentage of relative frequency')
    plt.xlabel('AA')
    plt.ylabel('%RF')
    plt.show()
    
    #Lets cluster by categorical properties using kmodes
    # One-hot encode the categorical data
    data_encoded = pd.get_dummies(df_final_data['Chemical Property'])
    
    # Number of clusters
    num_clusters = df_final_data['Chemical Property'].nunique()

    # Applying k-modes clustering
    km = KModes(n_clusters=num_clusters, init='Huang', n_init=5, verbose=1)
    clusters = km.fit_predict(data_encoded)

    # Adding cluster labels to your dataframe
    df_final_data['Cluster'] = clusters
    print(df_final_data)
    
    # Creating a crosstab for cluster and chemical properties
    cluster_crosstab = pd.crosstab(df_final_data['Cluster'], df_final_data['Chemical Property'])

    # Convert crosstab to a format suitable for Plotly
    cluster_crosstab = cluster_crosstab.reset_index().melt(id_vars='Cluster')

    # Creating an interactive stacked bar plot
    fig = px.bar(cluster_crosstab, x='Cluster', y='value', color='Chemical Property',
                title='Distribution of Chemical Properties in Each Cluster',
                labels={'value':'Count', 'Chemical Property':'Chemical Property'})
    fig.show()


    #Now lets save each cluster into a file for further analysis
    # Iterate over each cluster
    for cluster_number in df_final_data['Cluster'].unique():
        # Filter data for the current cluster
        cluster_data = df_final_data[df_final_data['Cluster'] == cluster_number]

        # Define the file name based on the cluster number
        file_name = f'{cluster_number}_members.csv'

        # Save the filtered data to a CSV file
        cluster_data.to_csv(file_name, index=False)

    



def main():
    # Check if help is requested
    if '-h' in sys.argv or '--help' in sys.argv:
        print_help()
        sys.exit(0)

    # Parsing command line arguments
    if len(sys.argv) != 7:
        print("Invalid arguments. Use -h or --help for usage instructions.")
        sys.exit(1)

    msa_file = ''
    aa_init = 0
    aa_pos = 0

    # Iterating over the arguments
    for i in range(1, len(sys.argv)):
        if sys.argv[i] == '-f' and i < len(sys.argv) - 1:
            msa_file = sys.argv[i + 1]
        elif sys.argv[i] == '-ai' and i < len(sys.argv) - 1:
            aa_init = int(sys.argv[i + 1])
        elif sys.argv[i] == '-ap' and i < len(sys.argv) - 1:
            aa_pos = int(sys.argv[i + 1])

    process(msa_file, aa_init, aa_pos)

def print_help():
    help_message = """
    Usage: python script.py -f [msafile] -ai [aa_init] -ap [aa_pos]

    Arguments:
    -f  [msafile] : Path to the multiple sequence alignment file. The file format should be as per the Anthony Nolan GitHub repository (https://github.com/ANHIG/IMGTHLA).
    -ai [aa_init] : Integer number indicating the initial amino acid position in the reference sequence.
    -ap [aa_pos]  : Integer specifying the position of the amino acids for analysis.

    Output:
    1. A plot of the relative percentage frequency of amino acid mutations.
    2. A plot showing the clustering of mutate AA based on their chemical properties.
    3. CSV files containing the data for each cluster.

    Description:
    This script processes a multiple sequence alignment file to analyze amino acid mutations. It categorizes amino acids based on their chemical properties and performs clustering. Outputs include visualizations and CSV files for further analysis.
    """
    print(help_message)

if __name__ == '__main__':
    main()
