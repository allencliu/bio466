def foo(adict, alist):
    my_d = {}
    for i, j in adict.items():
        print("key[{0}] value({1})".format(i, j))
        my_d[j] = i
    print("foo()", my_d)
    for i, j in enumerate(alist):
        adict[j]=adict[j]+i
    return my_d


dict1 = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
print("main 1 dict1 |{0}|".format(dict1))
dict2 = foo(dict1, ['A', 'T', 'G', 'C'])
print("main 2 dict1 |{0}|".format(dict1))
print("main 3 dict2 |{0}|".format(dict2))

myStr = 'ATGC'
for i in range(len(myStr)):
    print("For loop ({0}) |{3}| [{2}] <{1}>".format(i, myStr[i], dict1[myStr[i]], dict1[dict2[i]]))

num = 3
while num >= 1:
    num = num-1
    print("While loop ({1}) ({0})".format(myStr[num], myStr[num-1]))
    myStr = myStr+myStr[num]+myStr[num-1]

print("Finally,myStr=[{0}]".format(myStr))
