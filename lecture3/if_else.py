# Filename:if_else.py


def nucleotide_counter(X):
    """ Count nucleotide numbers for a given input sequence
    :para X:
    :return: no return
    """
    (A_count,C_count,G_count,T_count)=(0,0,0,0)
    if 'A' in X:
        A_count = X.count('A')
    elif 'G' in X:
        G_count = X.count('G')
    elif 'C' in X:
        C_count = X.count('C')
    elif 'T' in X:
        T_count = X.count('T')
    print("A={0}\nC={1}\nG={2}\nT={3}\n".format(A_count,C_count,G_count,T_count))


seq = input('Enter your sequence:\n')
nucleotide_counter(seq)
