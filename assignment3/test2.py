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
    print(aligned_seq1)
    print(aligned_seq2)
    
    
    # Initialize variables for match, mismatch, insertion, and deletion counts
    match_count = 0
    mismatch_count = 0
    insertion_count = 0
    deletion_count = 0

    alignment_line = ""  # Initialize the middle line of the alignment
    cigar_string = ""    # Initialize the CIGAR string

    # Calculate match, mismatch, insertion, and deletion counts
    for char1, char2 in zip(aligned_seq1, aligned_seq2):
        if char1 == char2:
            match_count += 1
            alignment_line += "|"  # Match symbol
            cigar_string += "M"
        elif char1 == '-' or char2 == '-':
            if char1 == '-':
                insertion_count += 1
                alignment_line += " "  # Insertion symbol
                cigar_string += "I"
            else:
                deletion_count += 1
                alignment_line += " "  # Deletion symbol
                cigar_string += "D"
        else:
            mismatch_count += 1
            alignment_line += " "  # Mismatch symbol
            cigar_string += "M"

    # Calculate match, mismatch, insertion, and deletion ratios
    total_length = len(aligned_seq1)
    match_ratio = match_count / total_length
    mismatch_ratio = mismatch_count / total_length
    insertion_ratio = insertion_count / total_length
    deletion_ratio = deletion_count / total_length
    print(alignment_line)

# Example usage:
sequence_names = ["Seq1", "Seq2", "Seq3"]
sequences = ["ACGT", "AGCT", "ACCG"]

alignment(sequence_names, sequences)