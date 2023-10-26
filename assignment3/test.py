
import re


def detect_homopolymer(sequence):
    homopolymers = []  # List to store detected homopolymers

    for nuc in ['A', 'T', 'G', 'C']:
        pattern = f"{nuc}{{3}}(?:{nuc}[^{nuc}]{{0,1}}){{4}}{nuc}{{3}}"
        for match in re.finditer(pattern, sequence):
            start = match.start()
            end = match.end() - 1  # Adjust the end position
            length = end - start + 1
            homopolymers.append(f'{nuc} {start} {end + 1} {length}')

    return homopolymers

print(detect_homopolymer("CGCCCAGCTGAAGTTCACCTTCCAGCCGCCTGAGGCGGCGATGCCGCTGAATGCAGCCCAGACAGACAGATTTTGGGGGCGGGCGCGCGCGCGCGCGGGGGCGGCGCGGGGCGCTTAAAATTGAAAACAAAAAACCATGC"))
print("CGCCCAGCTGAAGTTCACCTTCCAGCCGCCTGAGGCGGCGATGCCGCTGAATGCAGCCCAGACAGACAGATTTTGGGGGCGGGCGCGCGCGCGCGCGGGGGCGGCGCGGGGCGCTTAAAATTGAAAACAAAAAACCATGC"[97:110])


