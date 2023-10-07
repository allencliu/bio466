import liuac_a1 as bp

fasta = bp.read_fasta("Sequence.fasta")
seq = fasta[2]
print(seq[884:896])
print(len(seq[884:896]))