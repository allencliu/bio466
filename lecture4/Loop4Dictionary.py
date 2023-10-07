# Filename: Loops4Dictionary.py

dic = {'k1':10,'k2':11,'k3':12}
for key in dic.keys():
	print (key)
for val in dic.values():
	print(val)
for item in dic.items():
	print(item)

keyList=list(dic.keys())
print('keyList',keyList)
valueList=list(dic.values())
print('valueList',valueList)
itemList=list(dic.items())
print('itemList',itemList)

keySet=set(dic.keys())
print('keySet',keySet)
valueSet=set(dic.values())
print('valueSet',valueSet)
itemSet=set(dic.items())
print('itemSet',valueSet)
