# Filename: enumerate.py
# enumerate(iterable) will return both index and value

myStr = 'ATGC'
for i, j in enumerate(myStr):
    print(i, j)

myList = [1, 2, 3, 4, 5]
for x, y in enumerate(myList):
    print(x, y)

myTuple = (10, 20, 30, 40, 50)
for y, x in enumerate(myTuple):
    print(y, x)

mySet = {3, 3, 5, 5, 6, 6, 6, 7}
for m, n in enumerate(mySet):
    print(m, n)

myDict = {'A': 5, 'T': 3, 'G': 8, 'C': 10}
for j, i in enumerate(myDict):
    print(j, i, myDict[i])
