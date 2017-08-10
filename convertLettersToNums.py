openFile = open("genome.txt", "r")
targetFile = open("genomenums.txt","w")

for line in openFile.readlines():
    for char in line:
        if(char == "A"):
            targetFile.write("1\n")
        elif(char == "G"):
            targetFile.write("2\n")
        elif(char == "C"):
            targetFile.write("3\n")
        elif(char == "T"):
            targetFile.write("4\n")

openFile.close()
targetFile.close()
