X=''                 # X is an empty string
Line=0

with open('MySeq.fasta','r') as fh:
    next(fh)
    for n, line in enumerate(fh):
        Line-=n
        X+=line.strip()
        # print(X)
        # print(Line)
        
print(Line)