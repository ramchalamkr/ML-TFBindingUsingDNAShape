#To generate both the positive indices using DNAshapeR.

import itertools
from sklearn.feature_extraction.text import CountVectorizer
from sklearn import preprocessing
import numpy as np
import pandas as pd
import random
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import MultinomialNB
import pickle
import csv
from Bio import SeqIO
import subprocess


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

def generate_negative_kmer_indices(filePath):
    observationList = []
    count = 1
    #sequence = ''
    with open(filePath, 'r') as f:
        for i, line in enumerate(f):
            #sequence = sequence + line.rstrip()
            line = line.split()
            index1 = line[1]
            index2 = line[2]
            #if int(index1)-int(index2)>=40:
            observationList.append([line[0],int(index1),int(index1)+15])
    #observationList = observationList[1:]
    return observationList


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
                lowIndex = int(j[1])-20
                highIndex = int(j[2])+20
                fileName = 'chromFa/'+str(i)+'.fa'
                fastaSequence = read_fasta_file(fileName)
                testRun = fastaSequence[0][lowIndex:(highIndex)]
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
                
            
        #pickle.dump(mgwFileFinalList,open('ZNF263_mgw_Positive'+str(i)+'.p',"wb"))
        #pickle.dump(heltFileFinalList,open('ZNF263_helt_Positive'+str(i)+'.p',"wb"))
        #pickle.dump(protFileFinalList,open('ZNF263_prot_Positive'+str(i)+'.p',"wb"))
        #pickle.dump(rollFileFinalList,open('ZNF263_roll_Positive'+str(i)+'.p',"wb"))
        #pickle.dump(kmerList,open('UAK25_kmer_Positive'+str(i)+'.p',"wb"))
        #break
    return mgwFileFinalList,heltFileFinalList,protFileFinalList,rollFileFinalList,kmerList
#remove_comma_shape_file()


motifPosFileName = 'factorbookMotifPos.txt'
negativeMotifPosFileName = 'wgEncodeRegTfbsClusteredV3.GM12878.merged.bed'
targetTF = 'ZNF263'
kmerIndices = generate_kmer_indices(motifPosFileName,targetTF)
chrFileList = get_chrlist_TF(kmerIndices)

mgwPositiveList, heltPositiveList,protPositiveList,rollPositiveList,kmerList =  get_shape_list(targetTF,chrFileList,kmerIndices)
