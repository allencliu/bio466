
import regex


def detect_homopolymer(sequence):
    homopolymers = []  # List to store detected homopolymers

    for nuc in ['A', 'T', 'G', 'C']:
        # pattern = f"{nuc}{{3}}(?:{nuc}[^{nuc}]{{0,1}}){{4}}{nuc}{{3}}"
        
        pattern = r"({nuc}{3,}{nuc}{s<=1}{nuc}{3,})"
        pattern = pattern.replace("{nuc}", nuc)  # Replace the placeholder with the nucleotide

        regex_pattern = regex.compile(pattern)
        for match in regex_pattern.finditer(sequence):
            start = match.start()
            end = match.end() - 1  # Adjust the end position
            length = end - start + 1
            match_count = match.group(0).count(nuc)  # Count the occurrences of the nucleotide in the match
            if match_count >= 10:  # Check if the homopolymer length is at least 10
                match_count = match.group(0).count(nuc)  # Count the occurrences of the nucleotide in the match
                print(match_count)
                homopolymers.append(f'{nuc} {start} {end + 1} {length}')

    return homopolymers

print(detect_homopolymer("CGCCCAGCTGAAGTTCACCTTCCAGCCGCCTGAGGCGGCGATGCCGCTGAATGCAGCCCAGACAGACAGATTTTGGGGGCGGGCGCGCGCGCGCGCGGGGGCGGCGCGGGGCGCTTAAAATTGAAAACAAAAAACCATGC"))
print("CGCCCAGCTGAAGTTCACCTTCCAGCCGCCTGAGGCGGCGATGCCGCTGAATGCAGCCCAGACAGACAGATTTTGGGGGCGGGCGCGCGCGCGCGCGGGGGCGGCGCGGGGCGCTTAAAATTGAAAACAAAAAACCATGC"[123:134])


# res = regex.findall("(A{11,}){s<=1}", "CGCCCAGCTGAAGTTCACCTTCCAGCCGCCTGAGGCGGCGATGCCGCTGAATGCAGCCCAGACAGACAGATTTTGGGGGCGGGCGCGCGCGCGCGCGGGGGCGGCGCGGGGCGCTTAAAATTGAAAACAAAAAACCATGC", overlapped=True)
# res = regex.findall(r"(A{3,}A{s<=1}A{3,})", "CGCCCAGCTGAAGTTCACCTTCCAGCCGCCTGAGGCGGCGATGCCGCTGAATGCAGCCCAGACAGACAGATTTTGGGGGCGGGCGCGCGCGCGCGCGGGGGCGGCGCGGGGCGCTTAAAATTGAAAACAAAAAACCATGC", overlapped=True)
# print(res)