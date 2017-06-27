__author__ = 'djk'

import sys
f = open("IGVBatch.txt",'w')
sys.stdout = f

print("new\n"
      "genome\n"
      "load /Volumes/Seagate Backup Plus Drive/CCGD results_Nov2015/SCB0002h_TSC1-TSC2_B8/Bam files/igv_session_trial.xml\n"
      "snapshotDirectory Volumes/Seagate Backup Plus Drive/CCGD results_Nov2015/SCB0002h_TSC1-TSC2_B8/Bam files")
for line in open('calls.txt','r'):
    line2=line.replace("\n","")
    print("goto "+line2+"\n"
            "sort base\n"
            "maxPanelHeight 250\n"
            "snapshot "+line2+".png")
'''
print("load /Volumes/6TB-LaCie-7-12 1/EH07/EH07TN-16.xml\n"
      "snapshotDirectory /Volumes/6TB-LaCie-7-12 1/EH07/Snaps\n")
for line in open('calls.txt','r'):
    line2=line.replace("\n","")
    print("goto "+line2+"\n"
            "sort base\n"
            "maxPanelHeight 250\n"
            "snapshot "+line2+"-EH07T9-16.png")
'''