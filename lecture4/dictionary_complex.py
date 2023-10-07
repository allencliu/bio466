# Filename: dictionary_complex.py

MyD={}
MyD['This']=[0,1,2,3]
MyD['is']=('c','o','o','l')
MyD['Cool']={0:'c',1:'o',2:'o',3:'l'}

for item in MyD.items():
    print(item)

for key, value in MyD.items():
    print(key,value)
    if key == 'This' or key=='is':
        for x, y in enumerate(value):
            print("Index=[{0}] Value=[{1}]".format(x, y))
    else:
        for x, y in value.items():
            print(x, y, MyD['Cool'][x])
