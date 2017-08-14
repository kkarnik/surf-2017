# Import data from GenomAD Browser
source = open('exacalts.txt', 'r')
target = open('reformatexac.txt', 'w')

for line in source.readlines():
    if (line == 'A\n'):
        target.write('1\n')
    elif (line == 'G\n'):
        target.write('2\n')
    elif (line == 'C\n'):
        target.write('3\n')
    elif (line == 'T\n'):
        target.write('4\n')
    else:
        for char in line:
            if (char == 'A'):
                target.write('1')
            elif (char == 'G'):
                target.write('2')
            elif (char == 'C'):
                target.write('3')
            elif (char == 'T'):
                target.write('4')
        target.write('\n')

source.close()
target.close()