from sys import argv


def parse_config_file(filename):
    """Parses a configuration file and returns a dictionary of key-value pairs.

    Args:
        filename (str): The path/name to the configuration file to be parsed.

    Returns:
        dict: A dictionary containing key-value pairs from the configuration file (settings to be turned on or off).
    """
    config_dict = {}  # Initialize an empty dictionary to store configuration settings.
    fh = open(filename, "r")  # Open the configuration file in read mode.
    lines = fh.readlines()  # Read all lines from the file into a list.

    # Iterate through each line in the list of lines.
    for line in lines:
        parts = line.split(
            "="
        )  # Split each line into two parts using the "=" character as a delimiter.

        # Check if the line has exactly two parts (a key and a value).
        if len(parts) == 2:
            # Add the key-value pair to the config_dict dictionary, removing leading/trailing whitespace from the value.
            config_dict[parts[0]] = parts[1].strip()

    # Close the file.
    fh.close()

    return config_dict  # Return the dictionary containing the parsed configuration settings.


def read_fasta(filename):
    """Reads a FASTA file and extracts sequence information.

    Args:
        filename (str): The path/name of the FASTA file to be read.

    Returns:
        tuple: A tuple containing three lists: seq_name, seq_desc, and seq_base.
               - seq_name: List of sequence names.
               - seq_desc: List of sequence descriptions.
               - seq_base: List of sequence bases.
    """
    # Read the entire contents of the FASTA file into a string.
    fasta_file = open(filename, "r").read()

    # Split the file into a list of sequences using '>' as the delimiter.
    name_list = fasta_file.split(">")

    # Initialize empty lists to store sequence names, descriptions, and bases.
    seq_name = []
    seq_desc = []
    seq_base = []

    # Iterate through each sequence in the list (starting from index 1 to skip the empty first element).
    for seq in name_list[1:]:
        # Extract the sequence name from the first line.
        name = seq.split("\n")[0][:].split(" ")[0]

        # Extract the sequence description from the first line, excluding the name.
        description = " ".join(seq.split("\n")[0][1:].split(" ")[1:])

        # Extract the sequence bases by joining all lines after the first one.
        sequence = "".join(seq.split("\n")[1:])

        # Append the extracted information to their respective lists.
        seq_name.append(name)
        seq_desc.append(description)
        seq_base.append(sequence)

    # Return a tuple containing the three lists.
    return (seq_name, seq_desc, seq_base)


def print_in_fasta(seq_name, description, sequence, NucleotidesPerLine, spacer):
    """Prints sequence information in FASTA format.

    Args:
        seq_name (str): The name of the sequence.
        description (str): The description or additional information about the sequence.
        sequence (str): The sequence bases.
        NucleotidesPerLine (int): The maximum number of nucleotides to print per line.
        spacer (bool): A flag indicating whether to add spaces between nucleotides.
    """
    if (
        spacer == False
    ):  # If spacer is set to False, do not add spaces between nucleotides.
        print(">" + seq_name, description)  # Print the sequence name and description.
        for i in range(0, len(sequence), NucleotidesPerLine):
            # Print the sequence in lines of NucleotidesPerLine characters per line.
            print(sequence[i : i + NucleotidesPerLine])
    else:
        print(">" + seq_name, description)  # Print the sequence name and description.
        for i in range(0, len(sequence), NucleotidesPerLine):
            # Print the sequence with spaces added every 10 characters.
            chunk = sequence[i : i + NucleotidesPerLine]
            chunk_with_spaces = " ".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )
            print(chunk_with_spaces)


def print_with_ruler(seq_name, description, sequence, NucleotidesPerLine, spacer):
    """Prints sequence information with a ruler.

    Args:
        seq_name (str): The name of the sequence.
        description (str): The description or additional information about the sequence.
        sequence (str): The sequence bases.
        NucleotidesPerLine (int): The maximum number of nucleotides to print per line.
        spacer (bool): A flag indicating whether to add spaces between nucleotides.
    """
    if (
        spacer == False
    ):  # If spacer is set to False, do not add spaces between nucleotides.
        print(
            ">" + seq_name, description, "\n"
        )  # Print the sequence name and description.
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0) based on repeat_count.
        repeated_string = "".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        print(f"{1 :> 15}", end="")  # Print the first line number.
        for k in range(repeat_count - 1):
            print(f"{k + 2 :> 10}", end="")  # Print subsequent line numbers.
        print()  # Move to the next line.
        print(line_header)  # Print the ruler.

        # Print the sequence with line numbers and without spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Remove spaces between every 10 characters in the chunk.
            chunk_without_spaces = "".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_without_spaces)
    else:
        print(
            ">" + seq_name, description, "\n"
        )  # Print the sequence name and description.
        repeat_count = NucleotidesPerLine // 10

        # Generate a string of repeated digits from 1 to 10 (0 based on repeat_count.
        repeated_string = " ".join(["1234567890"] * repeat_count)

        # Create a line header with numbers and spaces.
        line_header = " ".join(["Line", repeated_string])

        # Print the header row with line numbers and ruler.
        print(f"{1 :> 15}", end="")  # Print the first line number.
        for k in range(repeat_count - 1):
            print(
                f"{k + 2 :> 11}", end=""
            )  # Print subsequent line numbers with extra spacing.
        print()  # Move to the next line.
        print(line_header)  # Print the ruler.

        # Print the sequence with line numbers and spaces.
        for i in range(0, len(sequence), NucleotidesPerLine):
            chunk = sequence[i : i + NucleotidesPerLine]

            # Add spaces between every 10 characters in the chunk.
            chunk_with_spaces = " ".join(
                [chunk[j : j + 10] for j in range(0, len(chunk), 10)]
            )

            # Print the line number and sequence chunk.
            print(f"{i // NucleotidesPerLine + 1 :> 4}", chunk_with_spaces)


def nucleotide_counter(sequence):
    """Count the occurrences of nucleotides in a DNA sequence.

    Args:
        sequence (str): The DNA sequence to be analyzed.

    Returns:
        tuple: A tuple containing the sequence length and counts of A, T, G, C, and N nucleotides.
    """
    # Get the length of the input DNA sequence.
    seq_length = len(sequence)

    # Initialize counters for each nucleotide type.
    A_count = 0
    T_count = 0
    G_count = 0
    C_count = 0
    N_count = 0

    # Iterate through each character in the sequence.
    for nucleotide in sequence:
        # Check the type of nucleotide and increment the corresponding counter.
        if nucleotide == "A":
            A_count += 1
        elif nucleotide == "T":
            T_count += 1
        elif nucleotide == "G":
            G_count += 1
        elif nucleotide == "C":
            C_count += 1
        elif nucleotide == "N":
            N_count += 1

    # Return a tuple containing the sequence length and nucleotide counts.
    return (seq_length, A_count, T_count, G_count, C_count, N_count)


def gc_content(sequence):
    """Calculate the GC content of a DNA sequence.

    Args:
        sequence (str): The DNA sequence for which GC content is calculated.

    Returns:
        float: The GC content as a fraction (between 0 and 1).
    """
    # Get the length of the input DNA sequence.
    seq_length = len(sequence)

    # Initialize counters for G and C nucleotides.
    G_count = 0
    C_count = 0

    # Iterate through each character in the sequence.
    for nucleotide in sequence:
        # Check if the nucleotide is 'G' (guanine) and increment the G count.
        if nucleotide == "G":
            G_count += 1
        # Check if the nucleotide is 'C' (cytosine) and increment the C count.
        elif nucleotide == "C":
            C_count += 1

    # Calculate the GC content as the sum of G and C counts divided by the sequence length.
    gc_content = (G_count + C_count) / seq_length

    # Return the calculated GC content as a fraction.
    return gc_content


def di_nucleotide_profile(sequence):
    """Calculate the counts of dinucleotides in a DNA sequence.

    Args:
        sequence (str): The DNA sequence for which dinucleotide counts are calculated.

    Returns:
        dict: A dictionary containing counts of dinucleotides as keys and their respective counts as values.
    """
    # Initialize a dictionary to store counts of dinucleotides.
    dinucleotide_counts = {
        "AA": 0,
        "AT": 0,
        "AG": 0,
        "AC": 0,
        "TT": 0,
        "TA": 0,
        "TC": 0,
        "TG": 0,
        "CC": 0,
        "CA": 0,
        "CG": 0,
        "CT": 0,
        "GG": 0,
        "GA": 0,
        "GT": 0,
        "GC": 0,
    }

    # Iterate through the sequence, considering pairs of consecutive nucleotides.
    for i in range(len(sequence) - 1):
        # Extract a dinucleotide pair from the sequence.
        dinucleotide = sequence[i : i + 2]

        # Check if the dinucleotide is present in the dinucleotide_counts dictionary.
        if dinucleotide in dinucleotide_counts:
            # Increment the count for the detected dinucleotide.
            dinucleotide_counts[dinucleotide] += 1

    # Return the dictionary containing counts of dinucleotides.
    return dinucleotide_counts


def cpg_island(sequence):
    """Identify CpG islands in a given DNA sequence.

    Args:
        sequence (str): The input DNA sequence.

    Returns:
        dict: A dictionary containing information about identified CpG islands.
            The keys are unique identifiers for each island (1-based index), and
            the values are strings in the format "start-end_length" where:
                - start: The start position of the CpG island.
                - end: The end position of the CpG island (exclusive).
                - length: The length of the CpG island.
    """
    i = 0
    count = 0
    cpg_dict = {}

    # Iterate through the sequence, considering pairs of consecutive nucleotides.
    while i in range(len(sequence) - 1):
        dinucleotide = sequence[i : i + 2]

        # Check if the current dinucleotide is "CG" (C followed by G).
        if dinucleotide == "CG":
            start = i  # Store the start position of the CpG island.

            # Continue checking for consecutive "CG" dinucleotides.
            while dinucleotide == "CG":
                dinucleotide = sequence[i : i + 2]
                i = i + 2  # Move to the next pair of nucleotides within the island.

            end = i - 2  # Store the end position of the CpG island.

            # Check if the island has a length of 6 or more.
            if end - start >= 6:
                count = count + 1

                # Create a range string (e.g., "start-end") and a length string.
                ranges = "-".join([str(start), str(end - 1)])
                value = "_".join([ranges, str(end - start)])

                # Store the CpG island information in the cpg_dict dictionary.
                cpg_dict[count] = value

        i = i + 2  # Move to the next pair of nucleotides in the sequence.

    # Return a dictionary containing information about identified CpG islands.
    return cpg_dict


def reverse_complement(sequence):
    """Calculate the reverse complementary DNA strand.

    Args:
        sequence (str): The input DNA sequence.

    Returns:
        str: The reverse complementary DNA strand.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_seq = sequence[::-1]
    reverse_comp_seq = "".join(complement[base] for base in reverse_seq)
    return reverse_comp_seq


def codon_profile(sequence):
    """Generate a codon profile for a given DNA sequence.

    Args:
        sequence (str): A string representing a DNA sequence.

    Returns:
        dict: A dictionary where keys are codons (3-letter DNA sequences) and values are initialized to 0.
    """
    bases = ["A", "T", "C", "G"]
    codon_dict = {}  # Initialize an empty dictionary to store the codon profile.
    for base1 in bases:
        for base2 in bases:
            for base3 in bases:
                codon = (
                    base1 + base2 + base3
                )  # Generate a codon by concatenating three bases.
                codon_dict[codon] = 0  # Initialize the codon count to 0.

    for i in range(len(sequence)):
        for codon in codon_dict:
            if sequence[i : i + 3] == codon:
                codon_dict[codon] = codon_dict[codon] + 1
    return codon_dict


def codon_profile_print(codon_dict):
    """Print the codon profile based on the provided codon dictionary.

    Args:
        codon_dict (dict): A dictionary containing codon counts, where the keys
            are codons (e.g., "TAA", "GCC") and the values are their respective counts.
    """
    print("Codon Profile:", end="\n\n")
    print(f"{'2nd':>20}")
    print(f"{'-'*35:>43}")
    print("1st\tT\tC\tA\tG\t3rd")
    # Define the order of bases
    bases = ["T", "C", "A", "G"]

    # Iterate over the first base
    for first_base in bases:
        row = first_base + "\t"

        # Counter to keep track of codons printed in the current line
        codons_printed = 0

        # Iterate over the second base
        for second_base in bases:
            # Iterate over the third base
            for third_base in bases:
                codon = first_base + third_base + second_base

                # Check if the codon exists in the dictionary
                if codon in codon_dict:
                    count = codon_dict[codon]

                # Append the codon and count to the row
                row += f"{codon}={count}\t"

                codons_printed += 1

                # Check if we have printed 4 codons, then start a new line
                if codons_printed == 4:
                    row += second_base  # Add the second base as 3rd column
                    print(row)
                    row = "\t"  # Start a new row
                    codons_printed = 0
        print()


def translation(dna_sequence):
    """Translate a DNA sequence into a protein sequence.

    Args:
        dna_sequence (str): The input DNA sequence to be translated.

    Returns:
        str: The resulting protein sequence.
    """
    trans_dic = {
        "UUU": "F",
        "UUC": "F",
        "UUA": "L",
        "UUG": "L",
        "CUU": "L",
        "CUC": "L",
        "CUA": "L",
        "CUG": "L",
        "AUU": "I",
        "AUC": "I",
        "AUA": "I",
        "AUG": "M",
        "GUU": "V",
        "GUC": "V",
        "GUA": "V",
        "GUG": "V",
        "UCU": "S",
        "UCC": "S",
        "UCA": "S",
        "UCG": "S",
        "CCU": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "ACU": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "GCU": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "UAU": "Y",
        "UAC": "Y",
        "UAA": "*",
        "UAG": "*",
        "CAU": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "AAU": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "GAU": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "UGU": "C",
        "UGC": "C",
        "UGA": "*",
        "UGG": "W",
        "CGU": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "AGU": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GGU": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }
    rna_sequence = dna_sequence.replace("T", "U")
    protein = []
    result = ""
    for i in range(0, len(rna_sequence) - 2, 3):
        protein.append(trans_dic[rna_sequence[i : i + 3]])
    result = result.join(protein)
    return result


def translation_6_frame(dna_sequence):
    """
    Translate a DNA sequence in all six reading frames.

    Args:
        dna_sequence (str): The input DNA sequence to be translated.

    Returns:
        Tuple of six protein sequences (str): The protein sequences translated from the six reading frames.
    """
    # Translate the original DNA sequence in the first reading frame
    protein_1 = translation(dna_sequence)

    # Translate the DNA sequence shifted by one base to the right (second reading frame)
    protein_2 = translation(dna_sequence[1:])

    # Translate the DNA sequence shifted by two bases to the right (third reading frame)
    protein_3 = translation(dna_sequence[2:])

    # Calculate the reverse complement of the DNA sequence
    rev_comp = reverse_complement(dna_sequence)
    # Translate the reverse complement in the first reading frame
    protein_4 = translation(rev_comp)

    # Translate the reverse complement shifted by one base to the right (fifth reading frame)
    protein_5 = translation(rev_comp[1:])

    # Translate the reverse complement shifted by two bases to the right (sixth reading frame)
    protein_6 = translation(rev_comp[2:])

    # Return a tuple containing the six protein sequences
    return (protein_1, protein_2, protein_3, protein_6[::-1], protein_5[::-1], protein_4[::-1])


def print_seq_fragment(seq_fragment, start, end):
    """
    Print a sequence fragment with annotations and translations.

    Args:
        seq_fragment (str): The sequence fragment to be printed.
        start (int): The start index of the fragment in the original sequence.
        end (int): The end index of the fragment in the original sequence.
    """
    
    # Calculate the length of the sequence fragment
    seq_length = end - start + 1

    # Extract the sequence from the fragment
    sequence = seq_fragment[start - 1 : end]

    # Translate the sequence in all six reading frames
    proteins = translation_6_frame(sequence)

    # Print the translations in the first three reading frames
    for p in proteins[3:]:
        print(format_seq_frag(sequence, p))

    # Print the reverse complement of the sequence
    print(reverse_complement(sequence)[::-1], end="\n")

    # Print a line of vertical bars to separate the annotations
    for i in range(seq_length):
        print("|", end="")
    print()

    # Print a line with sequence range annotation
    print(f"<{start}{'-' * (end - start - len(str(start)) - len(str(end)) - 1)}{end}>")

    # Print another line of vertical bars
    for i in range(seq_length):
        print("|", end="")
    print()

    # Print the original sequence
    print(sequence)
    
    # Print the translations in the last three reading frames
    for p in proteins[0:3]:
        print(format_seq_frag(sequence, p))
    print()
    
    # Perform additional analysis on the sequence (assumes analysis function exists)
    analysis(sequence)


def gene_analysis(gene_control, gene_target):
    """Analyze gene expression levels and determine their status.

    Args:
        gene_control (str): The expression level of the control gene.
        gene_target (str): The expression level of the target gene.

    Returns:
        str: A string representing the analysis result in the format "target:control:status",
            where status can be "*", "+", "-", or ".".
    """
    # Convert gene_control and gene_target to integers for comparison
    control_expression = int(gene_control)
    target_expression = int(gene_target)

    # Check if either gene expression is zero, indicating no expression
    if target_expression == 0 or control_expression == 0:
        return ":".join([gene_target, gene_control, "*"])

    # Calculate the ratio of target gene expression to control gene expression
    expression_ratio = target_expression / control_expression

    # Determine the status based on the expression ratio
    if expression_ratio >= 1.5:
        return ":".join([gene_target, gene_control, "+"])
    elif expression_ratio < (2 / 3):
        return ":".join([gene_target, gene_control, "-"])
    else:
        return ":".join([gene_target, gene_control, "."])


def get_expression(expression_file):
    """Read gene expression data from a file and create dictionaries for lung and prostate data.

    Args:
        expression_file (str): The path to the expression data file.

    Returns:
        tuple: A tuple containing two dictionaries:
            - dict_lung (dict): A dictionary mapping gene IDs to their expression analysis results in lung cancer.
            - dict_prostate (dict): A dictionary mapping gene IDs to their expression analysis results in prostate cancer.
    """
    # Initialize empty dictionaries for lung and prostate data
    dict_lung = {}
    dict_prostate = {}

    # Open the expression file for reading
    fh = open(expression_file, "r")

    # Read and discard the header line
    fh.readline()

    # Read the remaining lines in the file
    lines = fh.readlines()

    # Iterate over the lines and process gene expression data
    for line in lines:
        parts = line.split("\t")
        gene_id = parts[0]
        control_expression = parts[1]
        lung_expression = parts[2]
        prostate_expression = parts[3].strip()

        # Analyze gene expression for lung and prostate
        analysis_lung = gene_analysis(control_expression, lung_expression)
        analysis_prostate = gene_analysis(control_expression, prostate_expression)

        # Store the analysis results in the respective dictionaries
        dict_lung[gene_id] = analysis_lung
        dict_prostate[gene_id] = analysis_prostate

    # Close the file
    fh.close()

    # Return the dictionaries containing expression analysis results
    return dict_lung, dict_prostate


def print_gene_expression(dict_lung, dict_prostate):
    """Print differentially expressed genes detected for lung and prostate cancer-control comparisons.

    Args:
        dict_lung (dict): A dictionary mapping gene IDs to their expression analysis results in lung cancer.
        dict_prostate (dict): A dictionary mapping gene IDs to their expression analysis results in prostate cancer.
    """
    # Print differentially expressed genes for lung cancer-control comparison
    print(
        "1. Differentially expressed genes detected for lung cancer-control comparison:\n"
    )
    print("GeneID\tNumL:NumC:Exp")
    for key in dict_lung:
        print(f"{key} {dict_lung[key]}")
    print("\n")

    # Print differentially expressed genes for prostate cancer-control comparison
    print(
        "2. Differentially expressed genes detected for prostate cancer-control comparison:\n"
    )
    print("GeneID\tNumL:NumC:Exp")
    for key in dict_prostate:
        print(f"{key} {dict_prostate[key]}")
    print()


def gene_compare(dict_lung, dict_prostate):
    """Compare gene expression between lung and prostate cancer datasets.

    Args:
        dict_lung (dict): A dictionary mapping gene IDs to their expression analysis results in lung cancer.
        dict_prostate (dict): A dictionary mapping gene IDs to their expression analysis results in prostate cancer.

    Returns:
        tuple: A tuple containing four sorted lists:
            - common_genes (list): Genes expressed in both lung and prostate cancer.
            - lung_only_genes (list): Genes expressed only in lung cancer.
            - prostate_only_genes (list): Genes expressed only in prostate cancer.
            - neither_genes (list): Genes not expressed in both lung and prostate cancer.
    """
    # Extract sets of genes based on expression values
    genes_lung = {
        gene for gene, expression in dict_lung.items() if int(expression.split(":")[0]) > 0
    }
    genes_prostate = {
        gene
        for gene, expression in dict_prostate.items()
        if int(expression.split(":")[0]) > 0
    }

    # Extract sets of genes that end in '*' or '.'
    genes_lung_not_expressed = {
        gene for gene, expression in dict_lung.items() if int(expression.split(":")[0]) == 0
    }
    genes_prostate_not_expressed = {
        gene
        for gene, expression in dict_prostate.items()
        if int(expression.split(":")[0]) == 0
    }

    # Genes expressed in both lung and prostate cancer
    common_genes = sorted(list(genes_lung.intersection(genes_prostate)))

    # Genes expressed only in lung cancer
    lung_only_genes = sorted(list(genes_lung.difference(genes_prostate)))

    # Genes expressed only in prostate cancer
    prostate_only_genes = sorted(list(genes_prostate.difference(genes_lung)))

    # Genes expressed in neither lung nor prostate cancer (not expressed in both)
    neither_genes = sorted(
        list(genes_lung_not_expressed.intersection(genes_prostate_not_expressed))
    )

    return common_genes, lung_only_genes, prostate_only_genes, neither_genes


def print_gene_compare(
    common_genes, lung_only_genes, prostate_only_genes, neither_genes
):
    """Print four different groups of genes based on their expression patterns.

    Args:
        common_genes (list): Genes expressed in both lung and prostate cancer tissues.
        lung_only_genes (list): Genes expressed only in lung cancer tissues.
        prostate_only_genes (list): Genes expressed only in prostate cancer tissues.
        neither_genes (list): Genes expressed in neither lung nor prostate cancer tissues.
    """
    print("3. Four different groups of genes\n")

    # Print genes expressed in both lung and prostate cancer tissues
    print("[3.1] The genes expressed in both lung and prostate cancer tissues.\n")
    for gene in common_genes:
        print(gene, end=" ")
    print("\n")

    # Print genes expressed only in lung cancer tissues
    print("[3.2] The genes expressed only in lung cancer tissues.\n")
    for gene in lung_only_genes:
        print(gene, end=" ")
    print("\n")

    # Print genes expressed only in prostate cancer tissues
    print("[3.3] The genes expressed only in prostate cancer tissues.\n")
    for gene in prostate_only_genes:
        print(gene, end=" ")
    print("\n")

    # Print genes expressed in neither lung nor prostate cancer tissues
    print("[3.4] The genes expressed in neither lung nor prostate cancer tissues.\n")
    for gene in neither_genes:
        print(gene, end=" ")
    print()


def format_seq_frag(sequence, protein):
    """
    Distribute and center a protein within a sequence with padding.

    Args:
        sequence (str): The target string to center within.
        protein (str): The string to be centered.

    Returns:
        str: The centered string.
    """
    # Calculate the total number of spaces needed to center protein within sequence
    total_spaces = len(sequence) - len(protein) + 1

    # Calculate the number of spaces between characters to evenly distribute protein
    num_spaces_between_chars = total_spaces // (len(protein) - 1)

    # Initialize an empty string to store the distributed string
    distributed_protein = ""

    # Iterate through the characters of protein and add spaces between them
    for i, char in enumerate(protein):
        distributed_protein += char
        if i < len(protein) - 1:
            distributed_protein += " " * num_spaces_between_chars

    # Calculate the number of spaces needed on each side to center the distributed string
    total_padding_spaces = len(sequence) - len(distributed_protein)
    left_padding_spaces = total_padding_spaces // 2
    right_padding_spaces = total_padding_spaces - left_padding_spaces

    # Create the final centered string with padding
    centered_protein = (
        " " * left_padding_spaces + distributed_protein + " " * right_padding_spaces
    )

    return centered_protein


def analysis(sequence):
    """
    Perform various sequence analyses based on configuration settings.

    Args:
        sequence (str): The input sequence to be analyzed.
    """
    # Check if nucleotide counting is enabled in the configuration settings
    if config_dict["nucleotideCounter[Y|N]"] == "Y":
        # Calculate nucleotide counts
        n_tides = nucleotide_counter(sequence)
        # Print the nucleotide counts
        print(
            "Nucleotide Counts:",
            f"Seq Length={n_tides[0]}",
            f"A={n_tides[1]}",
            f"T={n_tides[2]}",
            f"G={n_tides[3]}",
            f"C={n_tides[4]}",
            f"N={n_tides[5]}",
        )

    # Check if GC content calculation is enabled in the configuration settings
    if config_dict["gcContent[Y|N]"] == "Y":
        # Calculate GC content and round it to two decimal places
        gc_cont = gc_content(sequence)
        # Print the GC content as a percentage
        print(f"GC Content={round(gc_cont * 100, 2)}%")

    # Check if dinucleotide profile calculation is enabled in the configuration settings
    if config_dict["dinucleotideProfil[Y|N]"] == "Y":
        # Calculate the dinucleotide profile
        di_profile = di_nucleotide_profile(sequence)
        # Print the dinucleotide profile
        print("Dinucleotide profile: ", end="")
        for dinucleotide in di_profile:
            print(f"{dinucleotide}={di_profile[dinucleotide]}", end=" ")
        print()

    # Check if CpG island prediction is enabled in the configuration settings
    if config_dict["CpGIsland[Y|N]"] == "Y":
        # Predict CpG islands
        cpg = cpg_island(sequence)
        # Print CpG islands
        print("CpG Islands: ", end="")
        for key in cpg:
            print("".join([str(key), "=", cpg[key]]), end=" ")
        print("\n")

    # Check if codon profile calculation is enabled in the configuration settings
    if "codonProfile[Y|N]" in config_dict and config_dict["codonProfile[Y|N]"] == "Y":
        # Calculate codon profile
        codon_profile_print(codon_profile(sequence))


def inquiry(config_dict, fasta_list):
    """Perform various analyses on selected sequences based on configuration settings.

    Args:
        config_dict (dict): A dictionary containing configuration settings.
        fasta_list (tuple): A tuple containing sequence data (names, descriptions, and sequences).

    """
    # Extract selected sequence indices and fragment start and end positions from config_dict.
    selected = [int(x) for x in config_dict["SelectedSeqs"].split(",")]
    start = int(config_dict["SeqFragmentStartPosition"])
    end = int(config_dict["SeqFragmentEndPostion"])
    seq_length = int(end) - int(start) + 1

    # Print information about the selected sequences and fragment positions.
    print(
        f"Among the {len(fasta_list[0])} sequences detected in the file: {config_dict['SeqFileName']}"
    )
    print(f"You have selected {selected} for the inquiry mode.")
    print(f"The start and end positions for sequence fragments: {start}-{end}")
    print()

    # Adjust selected sequence indices to be 0-based.
    selected = [x - 1 for x in selected]

    # Iterate through selected sequences.
    for seq in selected:
        # Check for invalid fragment positions.
        if start < 0 or end > len(fasta_list[2][seq]) or start > end:
            print(f"Invalid fragment for {fasta_list[0][seq]}!", end="\n\n")
        else:
            # Extract the sequence fragment.
            sequence = fasta_list[2][seq][start : end + 1]

            print(f">{fasta_list[0][seq]} {fasta_list[1][seq]}", end="\n\n")
            print(
                f"The selected fragment has a length of {seq_length} nucleotides:",
                end="\n\n",
            )

            if (
                "translation6Frames[Y|N]" in config_dict
                and config_dict["translation6Frames[Y|N]"] == "Y"
            ):
                print_seq_fragment(fasta_list[2][seq], start, end)
            else:
                # Print a representation of the sequence fragment with start and end markers.
                print(
                    f"<{start}{'-' * (end - start - len(str(start)) - len(str(end)) - 1)}{end}>"
                )
                for i in range(seq_length):
                    print("|", end="")
                print()

                # Print the sequence
                print(sequence, end="\n\n")
                analysis(sequence)


if __name__ == "__main__":
    """Main entry point for the script."""
    # Parse command line arguments and set configuration parameters.
    config_file = argv[1]  # Assuming command-line argument for the configuration file.
    config_dict = parse_config_file(config_file)
    print("Welcome Sequence Viewer!")
    print("Programmer:", config_dict["Programmer"])
    print("Email:", config_dict["Email"], "\n")

    # Read in the FASTA file containing sequence data.
    fasta_file = config_dict["SeqFileName"]
    fasta_list = read_fasta(fasta_file)

    # Display Mode
    if (
        config_dict["ViewSequenceInFastaFormat[Y|N]"] == "Y"
        and config_dict["DoYouNeedSpaceSeperator[Y|N]"] == "N"
    ):
        # Print sequences in FASTA format without spacing.
        for i in range(len(fasta_list[0])):
            group = fasta_list[0][i], fasta_list[1][i], fasta_list[2][i]
            print_in_fasta(
                *group, int(config_dict["NucleotidesPerLine[50|100]"]), False
            )
            print()
            # Analysis Mode
            analysis(fasta_list[2][i])

    elif (
        config_dict["ViewSequenceInFastaFormat[Y|N]"] == "Y"
        and config_dict["DoYouNeedSpaceSeperator[Y|N]"] == "Y"
    ):
        # Print sequences in FASTA format with spacing.
        for i in range(len(fasta_list[0])):
            group = fasta_list[0][i], fasta_list[1][i], fasta_list[2][i]
            print_in_fasta(*group, int(config_dict["NucleotidesPerLine[50|100]"]), True)
            print()
            analysis(fasta_list[2][i])

    elif (
        config_dict["ViewSequenceInFastaFormat[Y|N]"] == "N"
        and config_dict["DoYouNeedSpaceSeperator[Y|N]"] == "Y"
    ):
        # Print sequences with a ruler and spacing.
        for i in range(len(fasta_list[0])):
            group = fasta_list[0][i], fasta_list[1][i], fasta_list[2][i]
            print_with_ruler(
                *group, int(config_dict["NucleotidesPerLine[50|100]"]), True
            )
            print()
            analysis(fasta_list[2][i])

    elif (
        config_dict["ViewSequenceInFastaFormat[Y|N]"] == "N"
        and config_dict["DoYouNeedSpaceSeperator[Y|N]"] == "N"
    ):
        # Print sequences with a ruler without spacing.
        for i in range(len(fasta_list[0])):
            group = fasta_list[0][i], fasta_list[1][i], fasta_list[2][i]
            print_with_ruler(
                *group, int(config_dict["NucleotidesPerLine[50|100]"]), False
            )
            print()
            analysis(fasta_list[2][i])

    # Inquiry Mode
    inquiry(config_dict, fasta_list)

    dict_lung, dict_prostate = get_expression(config_dict["GeneExpFileName"])
    print_gene_expression(dict_lung, dict_prostate)
    common_genes, lung_only_genes, prostate_only_genes, neither_genes = gene_compare(
        dict_lung, dict_prostate
    )
    print_gene_compare(
        common_genes, lung_only_genes, prostate_only_genes, neither_genes
    )