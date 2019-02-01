# gets as input files from the R scirpt which have the cellID and the clusterID (files in input)
# prints the dataframe nicely to specified fileName *default is group_table

#outpu format: cellID-clusterID \n
#no whitespace

import numpy as np
import pandas as pd
from sklearn.preprocessing import scale

def readData(fileName,outfileName):
    fp = open( fileName, 'r')

    output_fileName = outfileName+"_group_table.txt"
    fout = open(output_fileName, "w") #clean the file
    fout = open(output_fileName, "a")

    line = fp.readline() #ignore headers
    cnt = 1
    array = []
    while line:
       line = fp.readline()
       if line == "":
           continue
       noq = 0 #number of quotation marks
       parsed_line = ""
       for c in line:
           if c== ' ':
               continue
           if c=='"':
               noq = noq+1
               continue
           if noq==3:
               parsed_line = parsed_line + c
       fout.write(parsed_line+'\n')
    fp.close()
    fout.close()


def main():
    pre_data = readData("files/preCellClusters.txt","patient9_pre");
    post_data = readData("files/postCellClusters.txt","patient9_post");

main()
