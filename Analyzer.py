__author__ = 'djk'

from Tkinter import *
import sys
import os
import subprocess
from multiprocessing.pool import ThreadPool as Pool
import shutil

#TSC1 9:135766735-135820020
#TSC2 16:2097990-2138713
#MTOR 1:11166588-11322608
#PTEN: 10:89623195-89728532
#PIK3CA: 3:178866311-178952497
#VHL: 3:10183319-10195354
#FLCN 17:17113000-17143000
#FH
class Application(Frame):
    # To run this application, the following files are needed:
    # numSamples.txt
    # genomewithinst.txt
    # merger.py
    # v12-q50.py
    # readcount.py
    # igv_session.xml
    # IGVBatch.py
    # Analyzer.py (which is this file itself)

    #@profile
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    #@profile
    def createWidgets(self):
        '''
        Creates the buttons and entry boxes for the user input for the analysis
        '''
        self.genbai = Button(self, text="Generate .bai", command=self.baigenerator)
        self.genbai.pack(side=TOP)

        L1 = Label(self, text="Location (chr:nt-nt):")
        L1.pack(side=LEFT)
        global E1
        E1 = Entry(self, bd =5)
        E1.pack(side=LEFT)

        L2 = Label(self, text="Minimum Variant Allele Frequency:")
        L2.pack(side=LEFT)
        global E2
        E2 = Entry(self, bd =5)
        E2.pack(side=LEFT)

        L3 = Label(self, text="Minimum read count in each direction:")
        L3.pack(side=LEFT)
        global E3
        E3 = Entry(self, bd =5)
        E3.pack(side=LEFT)

        self.start = Button(self, text="Begin Analysis", command=self.combinefunctions)
        self.start.pack()

        self.matlab = Button(self, text="Run matlab script", command=self.runmatlab)
        self.matlab.pack()


        #self.spreadsheet = Button(self, text="Excel Spreadsheet", command=self.createspreadsheet)

    def saveminfreq(self):
        '''
        Saves the input value for the minimum variant allele frequency into a text file, which will be imported as
        a variable in matlab.
        '''

        minfreq = E2.get()
        minfreqfile = open("minfreq.txt", "w")
        minfreqfile.write("%0.5f" % float(minfreq))
        minfreqfile.close()


    def saveminreadcount(self):
        '''
        Saves the input value for the minimum read count into a text file, which will be imported as a variable
        in matlab.
        '''

        minreadcount = E3.get()
        minreadcountfile = open("minreadcount.txt", "w")
        minreadcountfile.write("%d" % int(minreadcount))
        minreadcountfile.close()


    def runmatlab(self):
        '''
        Runs Matlab_commands_3_16_GenomADupdate.m script
        '''
        subprocess.call(['matlab -nojvm -r Matlab_commands_3_16_GenomADupdate'], shell=True)
        print('Matlab script completed. Open the file aafilterdata.csv\n')
        exit()

    #@profile
    def combinefunctions(self):
        '''
        Combines the functions that generate and modify the files in to a
        single function which can be invoked when the "Begin Analysis" button
        is pressed
        '''

        self.saveminfreq()
        self.saveminreadcount()

        self.namereader()

        if(not self.isvalidregion(geneloc)):
            print('Error, region is not valid\n')
            exit()
        else:
            print('Valid region.\n')

        self.pupgenerator()
        self.outgenerator()
        self.agenerator()
        self.zeroscreator()

        # Import data from GenomAD Browser
        source = open('genomeADalts.txt', 'r')
        target = open('reformatgenomeAD.txt', 'w')

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

        print("Python analysis complete. mergez.txt file and namelist.txt file generated.\n")

    def isvalidregion(self, geneloc):
        '''
        Compares the user input region values with the indices of the human genome reference values and makes sure
        that the start and end index are plus or minus 10000 of the reference
        '''
        i = geneloc.find(':') + 1
        j = geneloc.find('-')

        inputLen = len(geneloc)

        inputStart = int(geneloc[i : j])

        inputEnd = int(geneloc[j+1 : inputLen])

        refSeqFile = open("genomewithinst.txt", "r")

        fileData = refSeqFile.readline()

        indexbegin = fileData.find(':') + 1
        indexmid = fileData.find('-')
        indexend = fileData.find('\'') - 2

        refStart = int(fileData[indexbegin : indexmid])
        refEnd = int(fileData[indexmid+1 : indexend])

        text_file = open("initindex.txt", "w")
        text_file.write("%d" % refStart)
        text_file.close()

        target_file = open("genome.txt", "w")
        shutil.copyfileobj(refSeqFile, target_file)

        refSeqFile.close()
        target_file.close()

        return((refStart + 10000 == inputStart) and (refEnd - 10000 == inputEnd))

    #@profile
    def namereader(self):
        '''
        Reads through a list of files to get the bam files and stores the
        names of these files in namelist.txt
        '''
        global geneloc
        geneloc=E1.get()

        numSamples = 0

        oldstdout = sys.stdout
        f = open('namelist.txt','w')
        sys.stdout = f

        #prev file location was "/Users/Lana/Documents/Lana/BWH/Dr.Kwiatkowski/CCGD/Submission_Sep2015/Bam files
        #another old path was "/Users/guest/Desktop/SIS0009b"
        #path = "/Volumes/Untitled 1/TSC1-TSC2 SCB0002p"

        #path = "/Users/guest/Desktop/Karthik"

        path = os.path.realpath(__file__)

        path = path.replace('/Analyzer.py', '')

        dirs = os.listdir( path )

        for file in dirs:
            if file.endswith(".bam"):
                #file=file.replace(".pup","")
                print(file)
                numSamples += 1

        sys.stdout=oldstdout

        numSamplesFile = open("numSamples.txt", "w")
        numSamplesFile.write("%d" % numSamples)
        numSamplesFile.close()

    #@profile
    def baigenerator(self):
        '''
        Generates .bai files by calling a subprocess which calls the samtools
        index command on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['samtools index '+line3], shell=True)

    #@profile
    def pupgenerator(self):
        '''
        Generates .pup files by calling a subprocess which calls the samtools
        mpileup command on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['samtools mpileup -r '+geneloc+' '+line3+' > '+line3+'.pup'], shell=True)

    #@profile
    def outgenerator(self):
        '''
        Generates .out files by calling a subprocess on the UNIX command line
        that uses the Python script v12-q50.py
        '''

        # Previous version of the function was this:
        #for line in open('namelist.txt','r'):
        #    line2=line.replace("\n","")
        #    line3=line2.replace(" ",line2)
        #    subprocess.call(['python v12-q50.py '+line3+'.pup '+line3+'.out'], shell=True)

        lineList = open('namelist.txt', 'r')

        numSamplesFile = open("numSamples.txt", "r")
        numSamples = int(numSamplesFile.read())

        # Create the pool for the threads
        pool = Pool(numSamples)

        numSamplesFile.close()

        def out_loop_operation(line):
            '''
            Function that represents the operation done in the loop in the
            outgenerator function
            '''
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            #print("operation before call: %s\n" % line3)
            subprocess.call(['python v12-q50.py '+line3+'.pup '+line3+'.out'], shell=True)
            #print("operation after call: %s\n" % line3)



        # Make the threads run in parallel
        for line in lineList:
            #print("thread being run is %s\n" % line)
            pool.apply_async(out_loop_operation, (line,))

        pool.close()
        pool.join()



    #@profile
    def agenerator(self):
        '''
        Generates .a files by calling a subprocess which uses the cut command
        on the UNIX command line
        '''

        for line in open('namelist.txt','r'):
            line2 = line.replace("\n","")
            line3 = line2.replace(" ",line2)
            subprocess.call(['cut -f 2,3,6-9,11-14,16-18 '+line3+'.out > '+line3+'.a'], shell=True)

    #@profile

    #@profile
    def zeroscreator(self):
        '''
        Populates the mergez.txt file initially with zeroes
        '''
        #rowcount = str(subprocess.check_output(['wc -l < '+largestname], shell=True))

        subprocess.call(['python zeroscreator.py zeros.txt 70000 13'], shell=True)

        for line in open('namelist.txt','r'):
            line2 = line.replace("\n","")
            line3 = line2.replace(" ",line2)
            subprocess.call(['cat '+line3+'.a zeros.txt > '+line3+'.b'], shell=True)

        filenames=''
        for line in open('namelist.txt','r'):
            line2 = line.replace("\n","")
            line3 = line2.replace(" ",line2)
            filenames=filenames+line3+'.b '
        subprocess.call(['python merger.py '+filenames+'>>mergez.txt'], shell=True)

root = Tk()
app = Application(master=root)

app.master.title("Analyzer")

app.mainloop()
