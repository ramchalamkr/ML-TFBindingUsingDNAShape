# To generate the negatvie indices for all TFs for each of the chr files
import itertools
from sklearn.feature_extraction.text import CountVectorizer
from sklearn import preprocessing
import numpy as np
import pandas as pd
import random
import os
#from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import MultinomialNB
#from sklearn.neural_network import MLPClassifier
import pickle
import csv
from Bio import SeqIO
import subprocess
import re


def read_fasta_file(filePath):
    i = 0
    z = {}
    t = []

    for seq_record in SeqIO.parse(filePath, "fasta"):
        z[i] = seq_record.seq
        i+=1
    temp = str(z[0])
    t.append(temp)
    return t

def read_dnashape_file(filePath):
    i = 0
    z = {}
    from Bio import SeqIO
    for seq_record in SeqIO.parse(filePath, "fasta"):
        z[i] = seq_record.seq
        i+=1
    x = str(z[0])
    x1 = x.split(',')
    return x1[0:len(x1)-1]



def generate_kmer_indices(filePath,TFName):
    observationList = []
    count = 1
    #sequence = ''
    with open(filePath, 'r') as f:
        for i, line in enumerate(f):
            #sequence = sequence + line.rstrip()
            if(TFName in line):
                line = line.split()
                if(line[6] == '+'):
                    index1 = line[2]
                    index2 = line[3]
                    observationList.append([line[1],index1,index2])
    #observationList = observationList[1:]
    return observationList

def generate_negative_kmer_indices(filePath,chrList):
    observationList = []
    count = 1
    recordSize = 20
    #sequence = ''
    negativeKmerIndices = []
    for j in chrList:
        with open(filePath, 'r') as f:
            #print j
            for i, line in enumerate(f):
                line1=[]
                #print type(line)
                line1=re.split(r'\t+', line)
                if (line1[0]==j):
                    index1 = line1[1]
                    index2 = line1[2]
                    #print index2,index1
                    difference = int(index2)-int(index1)
                    noOfRecords = difference/recordSize
                    #print difference,noOfRecords
                    for i in range(noOfRecords):
                        lowerIndex = int(index1)
                        upperIndex = lowerIndex+20
                        #print lowerIndex,upperIndex
                        negativeKmerIndices.append([j,lowerIndex,upperIndex])
                        index1 = upperIndex
                    
    return negativeKmerIndices


def get_features(BindingIndices,files,l,l1,l2,l3,l4,flag):
    FinalBindingList =[]
    FinalBindingListShape =[]
    count = 0
    for [i,j,k] in BindingIndices:
        if i == files:
            count+=1
            FinalBindingList.append(l[0][int(j):int(k)].upper())
            FinalBindingListShape.append(l4[int(j):int(k)]) #,l2[int(j):int(k)],l3[int(j):int(k)],l4[int(j):int(k)]])
        if(count>=6800 and flag == 'N'):
            break
    return FinalBindingList ,FinalBindingListShape

def get_chrlist_TF(obvlist):
    UniqueList = []
    x = map(tuple,obvlist)
    seen = set()
    UniqueList = [item[0] for item in x if item[0] not in seen and not seen.add(item[0])]
    return UniqueList

def remove_comma_shape_file():
    chrFileList = ['chr9']
    for files in chrFileList:
        fileName = 'chromFa/'+files+'.fa'
        fileName2 = 'chromFa/'+files+'mod'+'.fa'
        with open(fileName+'.MGW') as infile, open(fileName2+'.MGW', 'w') as outfile:
            outfile.write(',\n'.join(infile.read().split('\n')) + '\n')
        with open(fileName+'.HelT') as infile, open(fileName2+'.HelT', 'w') as outfile:
            outfile.write(',\n'.join(infile.read().split('\n')) + '\n')
        with open(fileName+'.ProT') as infile, open(fileName2+'.ProT', 'w') as outfile:
            outfile.write(',\n'.join(infile.read().split('\n')) + '\n')
        with open(fileName+'.Roll') as infile, open(fileName2+'.Roll', 'w') as outfile:
            outfile.write(',\n'.join(infile.read().split('\n')) + '\n')

def add_comma_file(fileName):
    fileName2 = ''
    if 'MGW' in fileName:
        fileName2 = 'testMod.fa.MGW'
    elif 'HelT' in fileName:
        fileName2 = 'testMod.fa.HelT'
    elif 'ProT' in fileName:
        fileName2 = 'testMod.fa.ProT'
    elif 'Roll' in fileName:
        fileName2 = 'testMod.fa.Roll'
    #print fileName2
    with open(fileName) as infile, open(fileName2, 'w') as outfile:
            outfile.write(',\n'.join(infile.read().split('\n')) + '\n')
    return 0
#add_comma_file('test.fa.MGW')
def get_shape_list(tfName, chrList, indices ):
    for i in chrList:
        count = 0
        mgwFileFinalList = []
        heltFileFinalList = []
        protFileFinalList =[]
        rollFileFinalList =[]
        kmerList = []
        progressFile = open('progress.txt','w')
        progressFile.write(i)
        progressFile.close()
        for z,j in enumerate(indices):
            if i!=j[0]:
                continue
            #print i,j
            elif i==j[0]:
                #if (z%1==0):
                #    print 'processed',z,'records for',i
                lowIndex = int(j[1])-20
                highIndex = int(j[2])+20
                fileName = 'chromFa/'+str(i)+'.fa'
                fastaSequence = read_fasta_file(fileName)
                testRun = fastaSequence[0][lowIndex:(highIndex)]
                #print testRun
                #print fastaSequence[0][int(j[1]):int(j[2])]
                kmerList.append(testRun[20:35])
                #print len(testRun)
                fastaFileContent = '>seq'+"\n"+str(testRun)
                thefile = open('test.fa', 'w')
                thefile.write(fastaFileContent)
                thefile.close()
                s = ''
                s = "#!/usr/bin/Rscript" + "\n"
                s = s+ "library(DNAshapeR)" +"\n"
                s = s+ "getShape('test.fa')" +"\n"
                thefile = open('script.r', 'w')
                thefile.write(s)
                thefile.close()
                os.popen("./script.r").read()
                #(out,err) = proc.communicate()
                add_comma_file('test.fa.MGW')
                add_comma_file('test.fa.HelT')
                add_comma_file('test.fa.ProT')
                add_comma_file('test.fa.Roll')
                mgwFileList = read_dnashape_file('testMod.fa.MGW')
                heltFileList = read_dnashape_file('testMod.fa.HelT')
                protFileList = read_dnashape_file('testMod.fa.ProT')
                rollFileList = read_dnashape_file('testMod.fa.Roll')
                mgwFileList = mgwFileList[20:35]
                protFileList = protFileList[20:35]
                heltFileList = heltFileList[20:35]
                rollFileList = rollFileList[20:35]
                mgwFileFinalList.append(mgwFileList)
                heltFileFinalList.append(heltFileList)
                protFileFinalList.append(protFileList)
                rollFileFinalList.append(rollFileList)
                count+=1
                if count%300==0:
                    print 'Processed ',count,'records for',i
                
        pickle.dump(mgwFileFinalList,open('mgw_Negative'+str(i)+'.p',"wb"))
        pickle.dump(heltFileFinalList,open('helt_Negative'+str(i)+'.p',"wb"))
        pickle.dump(protFileFinalList,open('prot_Negative'+str(i)+'.p',"wb"))
        pickle.dump(rollFileFinalList,open('roll_Negative'+str(i)+'.p',"wb"))
        pickle.dump(kmerList,open('kmer_Negative'+str(i)+'.p',"wb"))
        break
    return mgwFileFinalList,heltFileFinalList,protFileFinalList,rollFileFinalList,kmerList



motifPosFileName = 'factorbookMotifPos.txt'
negativeMotifPosFileName = 'wgEncodeRegTfbsClusteredV3.GM12878.merged.bed'


chrFileList = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX']

negativeKmerIndices = generate_negative_kmer_indices(negativeMotifPosFileName,chrFileList)

pickle.dump(negativeKmerIndices,open('Negative_Indices.p',"wb"))
