from sys import argv

arg_num=len(argv)
print(arg_num)
print(type(argv))

for i in argv:
	print(i)

"""
try:
	f = open(argv[1],'r')
except IOError:
	print ("can't open {}".format(argv[1]))	
else:
	seq_list = []
	for num,line in enumerate(f):
		print ("[{0}]{1} [lineEnd]".format(num+1,line.strip()))
		if num+1 > 1:
			seq_list.append(line.strip())
	print ('Seq=',"".join(seq_list))
"""
# Names don't exist until bound. But the default int():0)