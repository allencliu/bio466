# Filename:loop_for_string.py


x = 'ATCG'   # x is an iterable object
for a in x:
    print(a, x)

for a, b in enumerate(x):
    print(a, b, x)

for a in range(len(x)):
    print(a,x[a],x)
