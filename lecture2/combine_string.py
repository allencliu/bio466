# Filename: combine_string.py

seq1 = 'AAATTT'
seq2 = 'CCCGGG'
combine = seq1 + seq2
combine = seq1 + ' ' + seq2
seq_size = len(combine)
print("{0} {1} {2} {3}".format(seq1,seq2,combine,seq_size))
print("({0})[{1}]|{2}|<{3}>".format(seq1,seq2,combine,seq_size))
