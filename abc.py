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
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

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

        #self.spreadsheet = Button(self, text="Excel Spreadsheet", command=self.createspreadsheet)

    def combinefunctions(self):
        '''
        Combines the functions that generate and modify the files in to a single function which can be invoked when
        the "Begin Analysis" button is pressed
        '''

        geneloc = self.namereader()
        self.pupgenerator(geneloc)
        self.outgenerator()
        self.agenerator()
        self.largestfile()
        self.zeroscreator()

    def namereader(self):
        '''
        Read through a list of files to get the bam files and store them in namelist.txt
        '''
        #global geneloc
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

    def baigenerator(self):
        '''
        Generates .bai files by calling a subprocess which calls the samtools index command on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['samtools index '+line3], shell=True)

    def pupgenerator(self, geneloc):
        '''
        Generates .pup files by calling a subprocess which calls the samtools mpileup command on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['samtools mpileup -r '+geneloc+' '+line3+' > '+line3+'.pup'], shell=True)

    def outgenerator(self):
        '''
        Generates .out files by calling a subprocess on the UNIX command line
        '''
        for line in open('namelist.txt','r'):
            line2=line.replace("\n","")
            line3=line2.replace(" ",line2)
            subprocess.call(['python v12-q50.py '+line3+'.pup '+line3+'.out'], shell=True)

    def agenerator(self):
        '''
        Generates .a files by calling a subprocess which uses the cut command on the UNIX command line
        '''

        for line in open('namelist.txt','r'):
            line2 = line.replace("\n","")
            line3 = line2.replace(" ",line2)
            subprocess.call(['cut -f 2,3,6-9,11-14,16-18 '+line3+'.out > '+line3+'.a'], shell=True)

    def largestfile(self):
        '''
        Finds the largest file in the set of .a files
        '''
        objects = os.listdir('.')

        sofar = 0
        global largestname
        largestname = ""

        for item in objects:
            if item.endswith('.a'):
                size = os.path.getsize(item)
                if size > sofar:
                        sofar = size
                        largestname = item

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

