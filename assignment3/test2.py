from NeedlemanWunsch_class_updated import Needleman_Wunsch
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

# # Example usage
# sequence = "ATGATGATGCATGGG"
# target = "ATG"
# motifs = motif_search(sequence, target)

# # Print the detected motifs
# for idx, motif_info in enumerate(motifs):
#     print(f"Motif {idx + 1}: {motif_info}")

# def modify_sequence(sequence, motifs):
#     modified_sequence = list(sequence)  # Convert the sequence to a list of characters for easy modification

#     for motif_info in motifs:
#         start_end, length = motif_info.split('_')
#         start, end = map(int, start_end.split('-'))

#         for i in range(start, end + 1):
#             modified_sequence[i] = sequence[i].lower()

#     return ''.join(modified_sequence)  # Convert the modified list back to a string

# # Example usage
# sequence = "ATGCTACATGCATATGCAGTCAATGCATACCCATGGG"
# motifs = motif_search(sequence, "ATG")

# modified_sequence = modify_sequence(sequence, motifs)
# print(modified_sequence)



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



# Example usage:
sequence_names = ["ENSEMBL008", "ENSEMBL009"]
sequences = ["AGGATATAGATAAGATAGACCACCCATAGA", "AGGAATAGATAAGATAGACAGGCGTATA"]

alignment(seq_name_list=sequence_names, seq_list=sequences)

