# Filename: range_index.py

myStr = 'ATGC'
for i in range(0, len(myStr)):
    print(i, myStr[i])

myList = [1, 2, 3, 4, 5]
for x in range(0, len(myList)):
    print(x, myList[x])

myTuple = (10, 20, 30, 40, 50)
for y in range(0, len(myTuple)):
    print(y, myTuple[y])

mySet = {3, 3, 5, 5, 6, 6, 6, 7}  # {3, 5, 6, 7}
myList = list(mySet)  # set does not support indexing
for m in range(0, len(myList)):
    print(m, myList[m])

myDict = {'A': 5, 'T': 3, 'G': 8, 'C': 10}  # dict does not support indexing
myKey = myDict.keys()
for i in myKey:
    print(i, myDict[i])
