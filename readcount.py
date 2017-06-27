__author__ = 'djk'

import sys
import os
import subprocess

for line in open('namelist.txt','r'):
        line2=line.replace("\n","")
        line3=line2.replace(" ",line2)
        for linea in open('exonlist.txt','r'):
                line4=linea.replace("\n","")
                line5=line4.replace(" ",line4)
                subprocess.call(['samtools view '+line3+' '+line5+' | wc -l'], shell=True)

