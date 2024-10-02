from Bio import SeqIO  # type: ignore

# Parse the fasta sequences from the file
fasta_sequences = SeqIO.parse("aligned.a2m", "fasta")

# List to store sequences that meet the criteria
updated_sequences = []

fit_CCCHCCHH = 0
others = 0

for fasta in fasta_sequences:
    # Remove lowercase characters from the sequence
    cleaned_seq = ''.join([char for char in str(fasta.seq) if not char.islower()])
    description_parts = fasta.description.split('|')
    original_description = '|'.join(description_parts[:4])  # Retain first 4 parts

    if (
        len(cleaned_seq) > 61
        and cleaned_seq[34] == 'C'
        and cleaned_seq[42] == 'C'
        and cleaned_seq[57] == 'C'
        and cleaned_seq[61] == 'H'
        and cleaned_seq[-33] == 'C'
        and cleaned_seq[-30] == 'C'
        and cleaned_seq[-14] == 'H'
        and cleaned_seq[-8] == 'H'
    ):
        # Sequence matches the CCCHCCHH pattern
        fit_CCCHCCHH += 1
        pattern = f'{cleaned_seq[34]}{cleaned_seq[42]}{cleaned_seq[57]}{cleaned_seq[61]}{cleaned_seq[-33]}{cleaned_seq[-30]}{cleaned_seq[-14]}{cleaned_seq[-8]}'
    else:
        # Sequence does not match
        others += 1
        pattern = f'{cleaned_seq[34]}{cleaned_seq[42]}{cleaned_seq[57]}{cleaned_seq[61]}{cleaned_seq[-33]}{cleaned_seq[-30]}{cleaned_seq[-14]}{cleaned_seq[-8]}'

    # Construct the new description
    new_description = f"{original_description}|{pattern}"

    # Update the description
    fasta.description = new_description

    updated_sequences.append(fasta)  # Add the updated fasta to the list

    # Optionally print the result for each sequence
    print(fasta.description)

# Summary
print(f"#fit_CCCHCCHH: {fit_CCCHCCHH}")
print(f"#others: {others}")

# Write updated sequences to a new fasta file manually
with open("updated_aligned.a2m", "w") as output_handle:
    for fasta in updated_sequences:
        # Write the header line
        output_handle.write(f">{fasta.description}\n")
        # Write the sequence
        output_handle.write(str(fasta.seq) + "\n")





# from Bio import SeqIO # type: ignore
# import pandas as pd # type: ignore

# # Parse the fasta sequences from the file
# fasta_sequences = SeqIO.parse("aligned.a2m", "fasta")

# # List to store sequences that meet the criteria
# updated_sequences = []

# fit_CCCHCCHH = 0
# others = 0

# for fasta in fasta_sequences:
#     # Remove lowercase characters from the sequence
#     cleaned_seq = ''.join([char for char in fasta.seq if not char.islower()])
#     description_parts = fasta.description.split('|')
#     original_description = '|'.join(description_parts[:4])

#     if len(cleaned_seq) > 61 and cleaned_seq[34] == 'C' and cleaned_seq[42] == 'C' and cleaned_seq[57] == 'C' and cleaned_seq[61] == 'H' and cleaned_seq[-33] == 'C' and cleaned_seq[-30] == 'C' and cleaned_seq[-14] == 'H' and cleaned_seq[-8] == 'H':
#         # Append "True" to the description of the fasta header
#         fit_CCCHCCHH += 1
#         new_description = f"{original_description}|{cleaned_seq[34]}{cleaned_seq[42]}{cleaned_seq[57]}{cleaned_seq[61]}{cleaned_seq[-33]}{cleaned_seq[-30]}{cleaned_seq[-14]}{cleaned_seq[-8]}"
#         # print(f"{fasta.description}: {cleaned_seq[34]},{cleaned_seq[42]},{cleaned_seq[57]},{cleaned_seq[61]},{cleaned_seq[-33]},{cleaned_seq[-30]},{cleaned_seq[-14]},{cleaned_seq[-8]} : {cleaned_seq}")
#     else:
#         others += 1
#         new_description += f"{original_description}|{cleaned_seq[34]}{cleaned_seq[42]}{cleaned_seq[57]}{cleaned_seq[61]}{cleaned_seq[-33]}{cleaned_seq[-30]}{cleaned_seq[-14]}{cleaned_seq[-8]}"

#     # print(fasta.description)
#     # new_description = fasta.description
#     # new_description += f"|{cleaned_seq[34]}{cleaned_seq[42]}{cleaned_seq[57]}{cleaned_seq[61]}{cleaned_seq[-33]}{cleaned_seq[-30]}{cleaned_seq[-14]}{cleaned_seq[-8]}"
#     fasta.description = new_description
#     print(fasta.description)
#     updated_sequences.append(fasta)  # Add the updated fasta to the list
    
# print(f"#fit_CCCHCCHH: {fit_CCCHCCHH}")
# print(f"#others: {others}")
# # Optionally, write updated sequences to a new fasta file
# with open("updated_aligned.a2m", "w") as output_handle:
#     SeqIO.write(updated_sequences, output_handle, "fasta")
