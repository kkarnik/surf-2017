# TODO: import the data from the genomAD browser and convert to format to be fed into matlab

source = open('genomeADalts.txt', 'r')
target = open('reformatgenomeAD.txt', 'w')

#print(source.readlines())

for line in source.readlines():
    if(line == 'A\n'):
        target.write('1\n')
    elif(line == 'G\n'):
        target.write('2\n')
    elif(line == 'C\n'):
        target.write('3\n')
    elif(line == 'T\n'):
        target.write('4\n')
    else:
        for char in line:
            if(char == 'A'):
                target.write('1')
            elif(char == 'G'):
                target.write('2')
            elif(char == 'C'):
                target.write('3')
            elif(char == 'T'):
                target.write('4')
        target.write('\n')

source.close()
target.close()
