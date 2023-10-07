# Filename:loop_while_string.py


x = 'AT'   # x is an iterable object
y = 'CG'
counter = 0
while counter <= 3:
    counter = counter + 1
    x = 'N' + x
    y = y + 'N'
    print("x={0} y={1} counter={2}".format(x,y,counter))

print("x={0} y={1} counter={2}".format(x,y,counter))