__author__ = 'djk'

from Tkinter import *
import sys
import os
import subprocess

#TSC1 9:135766735-135820020
#TSC2 16:2097990-2138713
#MTOR 1:11166588-11322608
#PTEN: 10:89623195-89728532
#PIK3CA: 3:178866311-178952497
#VHL: 3:10183319-10195354
#FLCN 17:17113000-17143000
#FH

class Application(Frame):
    @profile # For use with memory profiler (can be installed via pip)
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    @profile # For use with memory profiler (can be installed via pip)
    def createWidgets(self):
        self.genbai = Button(self, text="Generate .bai", command=self.baigenerator)
        self.genbai.pack(side=TOP)

        L1 = Label(self, text="Location (chr:nt-nt):")
        L1.pack(side=LEFT)
        global E1
        E1 = Entry(self, bd =5)
        E1.pack(side=LEFT)

        self.start = Button(self, text="Begin Analysis", command=self.combinefunctions)
        self.start.pack()

        self.matlab = Button(self, text="Run matlab script", command=self.runmatlab);
        self.matlab.pack();
        #self.spreadsheet = Button(self, text="Excel Spreadsheet", command=self.createspreadsheet)

    def runmatlab(self):
        '''
        Runs Matlab_commands_3_16_GenomADupdate.m script
        '''
        subprocess.call(['matlab -nojvm -r Matlab_commands_3_16_GenomADupdate'], shell=True)

    @profile # For use with memory profiler (can be installed via pip)
    def combinefunctions(self):
        '''
        Combines the functions that generate and modify the files in to a single function which can be invoked when
        the "Begin Analysis" button is pressed
        '''

        geneloc = self.namereader()
        self.pupgenerator(geneloc)
        self.outgenerator()
        self.agenerator()
        self.zeroscreator()

    @profile # For use with memory profiler (can be installed via pip)
    def namereader(self):
        '''
        Read through a list of files to get the bam files and store them in namelist.txt
        '''

        # Removed global variable in this file
        geneloc=E1.get()

        oldstdout = sys.stdout
        f = open('namelist.txt','w')
        sys.stdout = f

        #prev file location was "/Users/Lana/Documents/Lana/BWH/Dr.Kwiatkowski/CCGD/Submission_Sep2015/Bam files
        #another old path was "/Users/guest/Desktop/SIS0009b"
        path = "/Users/guest/Desktop/Karthik"

        dirs = os.listdir( path )

        for file in dirs:
            if file.endswith(".bam"):
                #file=file.replace(".pup","")
                print(file)

        sys.stdout=oldstdout

        return geneloc

    @profile # For use with memory profiler (can be installed via pip)
    def baigenerator(self):
        '''
        Generates .bai files by calling a subprocess which calls the samtools index command on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['samtools index '+line3], shell=True)

    @profile # For use with memory profiler (can be installed via pip)
    def pupgenerator(self, geneloc):
        '''
        Generates .pup files by calling a subprocess which calls the samtools mpileup command on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['samtools mpileup -r '+geneloc+' '+line3+' > '+line3+'.pup'], shell=True)

    @profile # For use with memory profiler (can be installed via pip)
    def outgenerator(self):
        '''
        Generates .out files by calling a subprocess on the UNIX command line that uses the Python script v12-q50.py
        '''

        # Previous version of the function was this:
        #for line in open('namelist.txt','r'):
        #    line2=line.replace("\n","")
        #    line3=line2.replace(" ",line2)
        #    subprocess.call(['python v12-q50.py '+line3+'.pup '+line3+'.out'], shell=True)

        lineList = open('namelist.txt', 'r')

        # Create the pool for the threads
        pool = Pool(2)

        def out_loop_operation(line):
            '''
            Function that represents the operation done in the loop in the outgenerator function
            '''
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['python v12-q50.py '+line3+'.pup '+line3+'.out'], shell=True)

        # Make the threads run in parallel
        for line in lineList:
            pool.apply_async(out_loop_operation, (line,))

        pool.close()
        pool.join()

    @profile # For use with memory profiler (can be installed via pip)
    def agenerator(self):
        '''
        Generates .a files by calling a subprocess which uses the cut command on the UNIX command line
        '''

        for line in open('namelist.txt','r'):
            line2 = line.replace("\n","")
            line3 = line2.replace(" ",line2)
            subprocess.call(['cut -f 2,3,6-9,11-14,16-18 '+line3+'.out > '+line3+'.a'], shell=True)

    @profile # For use with memory profiler (can be installed via pip)
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
