openFile = open("nt_coding_nums.txt", "r")
targetFile = open("codingLetters.txt","w")

for line in openFile.readlines():
    if(line == "1\n"):
        targetFile.write("A")
    elif(line == "2\n"):
        targetFile.write("G")
    elif(line == "3\n"):
        targetFile.write("C")
    elif(line == "4\n"):
        targetFile.write("T")

openFile.close()
targetFile.close()
