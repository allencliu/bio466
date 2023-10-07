# Filename:code_style.py
# This code shows two different coding styles
#      Function naming has different rules

def GetSum(X,Y):
    """ Calculate the sum for two numbers """
    sum=X+Y
    return sum

def get_sum(X,Y):
    """ Calculate the sum
        for two numbers
    """
    sum=X+Y
    return sum

result1=GetSum(2,3)
print(result1)

result2=get_sum(3,4)
print(result2)

print(help(GetSum))
print(help(get_sum))
