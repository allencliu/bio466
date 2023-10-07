import swalign
# choose your own values here… 2 and -1 are common.
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)

sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...
alignment = sw.align('ACACACTA','AGCACACA')
alignment.dump()

