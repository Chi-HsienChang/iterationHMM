import pandas as pd

# File paths
fasta_file_path = '../dataset/iteration2_to_iterationK/aligned_real.fasta'
tsv_file_path = '../dataset/iteration2_to_iterationK/LUC7p_species_phylogeny.tsv'

# Read the TSV file containing taxID and Named Lineage
tsv_df = pd.read_csv(tsv_file_path, sep='\t')
# Create a dictionary to map taxID to Named Lineage
taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}

# Initialize a dictionary to hold the mapping from FASTA name to Named Lineage
name_to_lineage = {}
name_to_taxID = {}
name_to_filterCCCHCCHH = {}

# Process the FASTA file
with open(fasta_file_path, 'r') as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            # Extract the identifier and remove the leading '>'
            identifier = line[1:]  # Remove '>' from the start of the header
            parts = identifier.split('|')
            name = parts[0]  # Assume the first part of the header is the name, now cleaned of '>'
            taxID = parts[-2].split(':')[-1]  # Extract taxID from the last part
            
            name_to_filterCCCHCCHH[name] = parts[-1]  # Extract taxID from the last part

            # Retrieve the Named Lineage using taxID
            named_lineage = taxid_to_lineage.get(taxID, 'NA,NA,NA,NA,NA')
            # Map the name to its corresponding Named Lineage
            name_to_lineage[name] = named_lineage
            print(named_lineage)
            name_to_taxID[name] = taxID


############################################################################################
############################################################################################
############################################################################################

fasta_file_path = 'i2_L3_after_cdhit.fasta'
output_file_path = 'updated_i2_L3_after_cdhit.fasta'



with open(fasta_file_path, 'r') as fasta_file, open(output_file_path, 'w') as output_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            parts = line[1:].split('|')  # Split the header into parts
            name = parts[0]  # Assume the first part is the name
            
            # Modify the header line to include name, lineage, and taxID
            line = f">{name}_{name_to_lineage[name].split(',')[3]}_{name_to_lineage[name].split(',')[4]}_{name_to_taxID[name]}"

        # Write the modified or unmodified line to the output file
        output_file.write(line + '\n')  # Add newline to preserve FASTA format

print("Updated FASTA file has been saved to:", output_file_path)
