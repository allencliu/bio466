# FileName: midterm-q2.py
def foo(seq):
    if len(seq)==1:
        result=seq
    else:
        result=seq[-1]+foo(seq[:-1])
    print("Inside foo()[{0}]".format(result))
    return result

condition = 1
seq = 'ATGC'
print("0 Before the loop: sequence has ({0})".format(seq))
seq = seq[0]+seq[-1]
num=len(seq)
print("1 Before the loop: sequence has <{0}>".format(seq))
while condition <= num:
    condition = condition+1
    seq = 'N'+seq
    seq1 = seq + 'N'
    print("condition=[{0}] seq=<{1}> seq1=({2})".format(condition, seq, seq1))

print("After loop condition=[{0}] seq=<{1}>".format(condition, seq))
seq=seq+'C'
foo(seq)

