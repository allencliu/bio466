# File Name: lessThanFilter.py

def lessThan(cutoffval,num_list):
    result=[]
    for num in num_list:
        if num < cutoffval:
            result.append(num)
    return result

def lessThanFilter(cutoffval,num_list):
    result1=[]
    result2=[]
    for num in num_list:
        if num < cutoffval:
            result1.append(num)
        else:
            result2.append(num)
    return result1,result2

print(lessThanFilter(10,[2,17,10,-3,42]))
print(lessThan(10,[2,17,10,-3,42]))