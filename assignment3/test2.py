from NeedlemanWunsch_class_updated import Needleman_Wunsch
import numpy as np
import matplotlib.pyplot as plt

def motif_search(sequence, target):
    detected_motifs = []
    seq_len = len(sequence)
    target_len = len(target)
    i = 0

    while i < seq_len:
        if sequence[i:i + target_len] == target:
            start_pos = i 
            end_pos = i + target_len - 1
            motif_length = target_len
            detected_motifs.append(f'{start_pos}-{end_pos}_{motif_length}')
            i += target_len
        else:
            i += 1

    return detected_motifs

def modify_sequence(sequence, motifs):
    modified_sequence = list(sequence)  # Convert the sequence to a list of characters for easy modification

    for motif_info in motifs:
        start_end, length = motif_info.split('_')
        start, end = map(int, start_end.split('-'))

        for i in range(start, end + 1):
            modified_sequence[i] = sequence[i].lower()

    return ''.join(modified_sequence)  # Convert the modified list back to a string

def alignment(seq_name_list, seq_list):
    # Ensure there are at least two sequences for alignment
    if len(seq_list) < 2:
        print("At least two sequences are required for alignment.")
        return

    # Define the NeedlemanWunsch object for sequence alignment
    # Align the first sequence with all other sequences
    first_sequence = seq_list[0]
    for i in range(1, len(seq_list)):
        aligner = Needleman_Wunsch(first_sequence, seq_list[i])
        aligned_seq1, aligned_seq2 = aligner.give_final_result()
        aligned_seq1 = seq_name_list[0] + "\t" + aligned_seq1
        aligned_seq2 = seq_name_list[i] + "\t" + aligned_seq2
        # Call the process_aligned_seq function to print alignments and CIGAR string
        process_aligned_seq(aligned_seq1, aligned_seq2)

def process_aligned_seq(aligned_seq1, aligned_seq2):
    def is_pyrimidine(char):
        return char in ('C', 'T')

    def is_purine(char):
        return char in ('A', 'G')
    
    print(aligned_seq1)
    # Initialize variables for match, mismatch, insertion, and deletion counts
    match_count = 0
    mismatch_count = 0
    insertion_count = 0
    deletion_count = 0
    pyr_pur_count = 0

    alignment_line = ""  # Initialize the middle line of the alignment
    cigar_string = ""    # Initialize the CIGAR string

    # Calculate match, mismatch, insertion, and deletion counts and build the CIGAR string
    current_operation = ""  # Track the current operation (M, I, D)
    current_count = 0  # Track the current count of contiguous operations
    
    for char1, char2 in zip(aligned_seq1.split("\t")[1], aligned_seq2.split("\t")[1]):
        if char1 == char2:
            match_count += 1
            alignment_line += "|"  # Match symbol
            if current_operation != "M":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "M"
                current_count = 1
            else:
                current_count += 1
        elif char1 == '-':
            insertion_count += 1
            alignment_line += " "  # Insertion symbol
            if current_operation != "I":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "I"
                current_count = 1
            else:
                current_count += 1
        elif char2 == '-':
            deletion_count += 1
            alignment_line += " "  # Deletion symbol
            if current_operation != "D":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "D"
                current_count = 1
            else:
                current_count += 1
        else:
            if is_pyrimidine(char1) and is_pyrimidine(char2):
                alignment_line += ":"  # Colon for pyrimidines
                pyr_pur_count += 1
            elif is_purine(char1) and is_purine(char2):
                alignment_line += ":"  # Colon for purines
                pyr_pur_count += 1
            else:
                alignment_line += " "  # Mismatch symbol
            mismatch_count += 1
            if current_operation != "M":
                if current_operation:
                    cigar_string += str(current_count) + current_operation
                current_operation = "M"
                current_count = 1
            else:
                current_count += 1

    # Append the last operation to the CIGAR string
    if current_operation:
        cigar_string += str(current_count) + current_operation

    # Calculate match, mismatch, insertion, and deletion ratios
    total_length = len(aligned_seq1)
    match_ratio = match_count / total_length
    mismatch_ratio = mismatch_count / total_length
    insertion_ratio = insertion_count / total_length
    deletion_ratio = deletion_count / total_length

    alignment_spaces = len(aligned_seq1.split("\t")[0]) * ' '
    # print(len(alignment_spaces))
    alignment_line = alignment_spaces + "\t" + alignment_line
    print(alignment_line)
    print(aligned_seq2)
    print(f"Identities={match_count}/{total_length} ({match_ratio * 100:.1f}%)", end= " ")
    print(f"Similarity={(match_count + pyr_pur_count)}/{total_length} ({(match_count + pyr_pur_count) / total_length * 100:.1f}%)", end=" ")
    print(f"Gaps={insertion_count + deletion_count}/{total_length} ({(insertion_count + deletion_count) / total_length * 100:.1f}%)", end=" ")
    print(f"CIGAR = {cigar_string}", end="\n")

def positional_matrix(seq_list):
    # Get the length of the sequences and the alphabet size
    seq_length = len(seq_list[0])
    nucleotides = ['A', 'T', 'G', 'C']
    nucleotides_size = len(nucleotides)

    # Initialize an empty positional matrix filled with zeros
    pos_matrix = np.zeros((nucleotides_size, seq_length), dtype=int)

    # Create a mapping from characters to column indices in the matrix
    char_to_index = {'A': 0, 'T': 1, 'G': 2, 'C': 3}

    # Iterate through the sequences character by character
    for sequence in seq_list:
        for position, char in enumerate(sequence):
            if char in char_to_index:
                row_index = char_to_index[char]
                pos_matrix[row_index, position] += 1

    return pos_matrix

def show_positional_matrix(my_array):
    # Get the number of positions and alphabet size from the positional matrix
    nucleotide_size, num_positions = my_array.shape
    left_align = len(str(num_positions))

    # Print position headers
    print("Position", end=" ")
    for i in range(1, num_positions + 1):
        print("{: <{width}}".format(i, width=left_align), end=" ")
    print()

    # Print a line separator
    print("--------", end=" ")
    for i in range(num_positions):
        print(f"{left_align * '-'}", end=" ")
    print()

    # Print the character counts for each nucleotide
    nucleotides = ['A', 'T', 'G', 'C']
    for row_index, nucleotide in enumerate(nucleotides):
        print(nucleotide, end=f"{8 * ' '}")
        for position in range(num_positions):
            print(f"{my_array[row_index, position]}", end=f"{left_align * ' '}")
        print()
        
        
    
    # Create a color map for nucleotides
    colors = {'A': 'blue', 'T': 'red', 'G': 'gray', 'C': 'orange'}
    
    # Calculate nucleotide percentages
    total_counts = np.sum(my_array, axis=0)
    percentages = (my_array / total_counts) * 100

    plt.figure(figsize=(16, 9))
    # Plot a stacked bar graph
    positions = range(1, num_positions + 1)
    nucleotides = ['A', 'T', 'G', 'C']
    bottom = np.zeros(num_positions)  # Initialize the bottom positions for stacking

    for nucleotide in nucleotides:
        counts = percentages[nucleotides.index(nucleotide)]
        plt.bar(positions, counts, label=nucleotide, color=colors[nucleotide], bottom=bottom)
        bottom += counts

    plt.title('Positional Nucleotide Matrix')
    
    plt.yticks(range(0, 101, 10))
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False, borderaxespad=0.5)
    # Save the stacked bar chart as a PNG file
    plt.savefig('PNM.jpg', format='JPG')

if __name__ == "__main__":
    
    # Example usage:
    # sequence_names = ["ENSEMBL008", "ENSEMBL009"]
    # sequences = ["AGGATATAGATAAGATAGACCACCCATAGA", "AGGAATAGATAAGATAGACAGGCGTATA"]

    # alignment(seq_name_list=sequence_names, seq_list=sequences)
    # Example usage:
    sequences = ["TGTCCAACGGGCCGAGGTTGTCTCTTTCGAGATCTTGTCGCGGGGGGGGGCTGCCTGTG",
                 "AGTGCGGGTCACAGATGCCTCGCACCCCTCCCCCGACAATGTGGCCCGTATGGAGGGTC",
                 "CCTGCCGCCGCGGGCGCGCATGCGCCCCGGCCCGTACTGGCCGGTGGTGCACCCGGCTG",
                 "CCGACACCCCGAGCGGGCCCGGGTTTTCACGTGCCTGGTCCCGCCACGCCCACCACATC"]
    result = positional_matrix(sequences)
    # print(result)
    show_positional_matrix(result)