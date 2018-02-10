#To generate the positive features for MAX TF
import itertools 
from sklearn import preprocessing
import numpy as np
import random
import pickle
import re
import math
import pickle
from sklearn.cross_validation import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_predict
from sklearn import metrics
from sklearn.ensemble import ExtraTreesClassifier

def generate_feature_set(shapelist):
	x = np.zeros(20)
	#for i in shapelist:
	for j in shapelist:
		x = np.vstack([x,j])
	x = x[1:]
	return x


def get_pwm_score(shapelist,M):
	x = np.zeros(1)
	for i in shapelist:
		for z in i:
			p = 1
			for j,k in enumerate(z):
				if(k=='A'):
					p*=M[0,j]
				elif(k=='C'):
					p*=M[1,j]
				elif(k=='G'):
					p*=M[2,j]
				elif(k=='T'):
					p*=M[3,j]
				else:
					p*=0
			x = np.vstack([x,p])
	x = x[1:]
	return x

def generate_feature_set_kmer(shapelist):
	x = np.zeros(20)
	for i in shapelist:
		#print i
		for z in i:
			#print z
			p = []
			for j,k in enumerate(z):
				if(k=='A' or k == 'a'):
					p.append("1000")
				elif(k=='C' or k == 'c'):
					p.append("0100")
				elif(k=='G' or k == 'g'):
					p.append("0010")
				elif(k=='T' or k == 't'):
					p.append("0001")
				else:
					p.append("abcd")
			x = np.vstack([x,p])
	x = x[1:]
	return x

a1 = pickle.load(open('featureset/chr1positiveexamples/heltchr1only20.p',"rb"))
a2 = pickle.load(open('featureset/chr2positivefiles/heltchr2only20.p',"rb"))
a3 = pickle.load(open('featureset/chr3positivefiles/heltchr3only20.p',"rb"))
a4 = pickle.load(open('featureset/chr4positivefiles/heltchr4only20.p',"rb"))
a5 = pickle.load(open('featureset/chr5positivefiles/heltchr5only20.p',"rb"))
a6 = pickle.load(open('featureset/chr6positivefiles/heltchr6only20.p',"rb"))
a7 = pickle.load(open('featureset/chr7positivefiles/heltchr7only20.p',"rb"))
a8 = pickle.load(open('featureset/chr8positivefiles/heltchr8only20.p',"rb"))
a9 = pickle.load(open('featureset/shapeFeatureFinal10helt.p',"rb"))
a10 = pickle.load(open('featureset/shapeFeatureFinalhelt9.p',"rb"))
a11= pickle.load(open('featureset/shapeFeatureFinalhelt11.p',"rb"))
a12 = pickle.load(open('featureset/shapeFeatureFinalhelt12.p',"rb"))
a13= pickle.load(open('featureset/shapeFeatureFinalhelt13.p',"rb"))
a14= pickle.load(open('featureset/shapeFeatureFinalhelt14.p',"rb"))
a15= pickle.load(open('featureset/shapeFeatureFinalhelt1518.p',"rb"))
a16= pickle.load(open('featureset/shapeFeatureFinalhelt19X.p',"rb"))

print 'loaded helt pickle files'

b1 = pickle.load(open('featureset/chr1positiveexamples/mgwchr1only20.p',"rb"))
b2 = pickle.load(open('featureset/chr2positivefiles/mgwchr2only20.p',"rb"))
b3 = pickle.load(open('featureset/chr3positivefiles/mgwchr3only20.p',"rb"))
b4 = pickle.load(open('featureset/chr4positivefiles/mgwchr4only20.p',"rb"))
b5 = pickle.load(open('featureset/chr5positivefiles/mgwchr5only20.p',"rb"))
b6 = pickle.load(open('featureset/chr6positivefiles/mgwchr6only20.p',"rb"))
b7 = pickle.load(open('featureset/chr7positivefiles/mgwchr7only20.p',"rb"))
b8 = pickle.load(open('featureset/chr8positivefiles/mgwchr8only20.p',"rb"))
b9 = pickle.load(open('featureset/shapeFeatureFinal10.p',"rb"))
b10 = pickle.load(open('featureset/shapeFeatureFinal9.p',"rb"))
b11= pickle.load(open('featureset/shapeFeatureFinal11.p',"rb"))
b12 = pickle.load(open('featureset/shapeFeatureFinal12.p',"rb"))
b13= pickle.load(open('featureset/shapeFeatureFinal13.p',"rb"))
b14= pickle.load(open('featureset/shapeFeatureFinal14.p',"rb"))
b15= pickle.load(open('featureset/shapeFeatureFinal1518.p',"rb"))
b16= pickle.load(open('featureset/shapeFeatureFinal19X.p',"rb"))

print 'loaded mgw pickle files'

c1 = pickle.load(open('featureset/chr1positiveexamples/protchr1only20.p',"rb"))
c2 = pickle.load(open('featureset/chr2positivefiles/protchr2only20.p',"rb"))
c3 = pickle.load(open('featureset/chr3positivefiles/protchr3only20.p',"rb"))
c4 = pickle.load(open('featureset/chr4positivefiles/protchr4only20.p',"rb"))
c5 = pickle.load(open('featureset/chr5positivefiles/protchr5only20.p',"rb"))
c6 = pickle.load(open('featureset/chr6positivefiles/protchr6only20.p',"rb"))
c7 = pickle.load(open('featureset/chr7positivefiles/protchr7only20.p',"rb"))
c8 = pickle.load(open('featureset/chr8positivefiles/protchr8only20.p',"rb"))
c9 = pickle.load(open('featureset/shapeFeatureFinal10prot.p',"rb"))
c10 = pickle.load(open('featureset/shapeFeatureFinalprot9.p',"rb"))
c11= pickle.load(open('featureset/shapeFeatureFinalprot11.p',"rb"))
c12 = pickle.load(open('featureset/shapeFeatureFinalprot12.p',"rb"))
c13= pickle.load(open('featureset/shapeFeatureFinalprot13.p',"rb"))
c14= pickle.load(open('featureset/shapeFeatureFinalprot14.p',"rb"))
c15= pickle.load(open('featureset/shapeFeatureFinalprot1518.p',"rb"))
c16= pickle.load(open('featureset/shapeFeatureFinalprot19X.p',"rb"))

print 'loaded prot pickle files'

d1 = pickle.load(open('featureset/chr1positiveexamples/rollchr1only20.p',"rb"))
d2 = pickle.load(open('featureset/chr2positivefiles/rollchr2only20.p',"rb"))
d3 = pickle.load(open('featureset/chr3positivefiles/rollchr3only20.p',"rb"))
d4 = pickle.load(open('featureset/chr4positivefiles/rollchr4only20.p',"rb"))
d5 = pickle.load(open('featureset/chr5positivefiles/rollchr5only20.p',"rb"))
d6 = pickle.load(open('featureset/chr6positivefiles/rollchr6only20.p',"rb"))
d7 = pickle.load(open('featureset/chr7positivefiles/rollchr7only20.p',"rb"))
d8 = pickle.load(open('featureset/chr8positivefiles/rollchr8only20.p',"rb"))
d9 = pickle.load(open('featureset/shapeFeatureFinalroll10.p',"rb"))
d10 = pickle.load(open('featureset/shapeFeatureFinalroll9.p',"rb"))
d11= pickle.load(open('featureset/shapeFeatureFinalroll11.p',"rb"))
d12 = pickle.load(open('featureset/shapeFeatureFinalroll12.p',"rb"))
d13= pickle.load(open('featureset/shapeFeatureFinalroll13.p',"rb"))
d14= pickle.load(open('featureset/shapeFeatureFinalroll14.p',"rb"))
d15= pickle.load(open('featureset/shapeFeatureFinalroll1518.p',"rb"))
d16= pickle.load(open('featureset/shapeFeatureFinalroll19X.p',"rb"))

print 'loaded roll pickle files'

e1Old = pickle.load(open('featureset/chr1positiveexamples/kmerchr1only20.p',"rb"))
e1 =[]
e1.append(e1Old)
#print e1
e2Old = pickle.load(open('featureset/chr2positivefiles/kmerchr2only20.p',"rb"))
e2 =[]
e2.append(e2Old)
e3Old = pickle.load(open('featureset/chr3positivefiles/kmerchr3only20.p',"rb"))
e3 =[]
e3.append(e3Old)
e4Old = pickle.load(open('featureset/chr4positivefiles/kmerchr4only20.p',"rb"))
e4 =[]
e4.append(e4Old)
e5Old = pickle.load(open('featureset/chr5positivefiles/kmerchr5only20.p',"rb"))
e5 =[]
e5.append(e5Old)
e6Old = pickle.load(open('featureset/chr6positivefiles/kmerchr6only20.p',"rb"))
e6 =[]
e6.append(e6Old)
e7Old = pickle.load(open('featureset/chr7positivefiles/kmerchr7only20.p',"rb"))
e7 =[]
e7.append(e7Old)
e8Old = pickle.load(open('featureset/chr8positivefiles/kmerchr8only20.p',"rb"))
e8 =[]
e8.append(e8Old)
e9 = pickle.load(open('featureset/kmerFeatureFinal10.p',"rb"))
e10 = pickle.load(open('featureset/kmerFeatureFinal9.p',"rb"))
e11= pickle.load(open('featureset/kmerFeatureFinal11.p',"rb"))
e12 = pickle.load(open('featureset/kmerFeatureFinal12.p',"rb"))
e13= pickle.load(open('featureset/kmerFeatureFinal13.p',"rb"))
e14= pickle.load(open('featureset/kmerFeatureFinal14.p',"rb"))
e15= pickle.load(open('featureset/kmerFeatureFinal1518.p',"rb"))

e16= pickle.load(open('featureset/kmerFeatureFinal19X.p',"rb"))

print 'loaded kmer pickle files'
#print a1

x11 = generate_feature_set(a1)
#print x11
x12 = generate_feature_set(a2)
x13 = generate_feature_set(a3)
x14 = generate_feature_set(a4)
x15 = generate_feature_set(a5)
x16 = generate_feature_set(a6)
x17 = generate_feature_set(a7)
x18 = generate_feature_set(a8)
x19 = generate_feature_set(a9)
x110 = generate_feature_set(a10)

x111 = generate_feature_set(a11)
x112 = generate_feature_set(a12)
x113 = generate_feature_set(a13)
x114 = generate_feature_set(a14)
x115 = generate_feature_set(a15)
x116 = generate_feature_set(a16)


print " Helt read done"

x21 = generate_feature_set(b1)
x22 = generate_feature_set(b2)
x23 = generate_feature_set(b3)
x24 = generate_feature_set(b4)
x25 = generate_feature_set(b5)
x26 = generate_feature_set(b6)
x27 = generate_feature_set(b7)
x28 = generate_feature_set(b8)
x29 = generate_feature_set(b9)
x210 = generate_feature_set(b10)
x211 = generate_feature_set(b11)
x212 = generate_feature_set(b12)
x213 = generate_feature_set(b13)
x214 = generate_feature_set(b14)
x215 = generate_feature_set(b15)
x216 = generate_feature_set(b16)

print " MGW read done"


x31 = generate_feature_set(c1)
x32 = generate_feature_set(c2)
x33 = generate_feature_set(c3)
x34 = generate_feature_set(c4)
x35 = generate_feature_set(c5)
x36 = generate_feature_set(c6)
x37 = generate_feature_set(c7)
x38 = generate_feature_set(c8)
x39 = generate_feature_set(c9)
x310 = generate_feature_set(c10)
x311 = generate_feature_set(c11)
x312 = generate_feature_set(c12)
x313 = generate_feature_set(c13)
x314 = generate_feature_set(c14)
x315 = generate_feature_set(c15)
x316 = generate_feature_set(c16)


print " prot read done"

x41 = generate_feature_set(d1)
x42 = generate_feature_set(d2)
x43 = generate_feature_set(d3)
x44 = generate_feature_set(d4)
x45 = generate_feature_set(d5)
x46 = generate_feature_set(d6)
x47 = generate_feature_set(d7)
x48 = generate_feature_set(d8)
x49 = generate_feature_set(d9)
x410 = generate_feature_set(d10)
x411 = generate_feature_set(d11)
x412 = generate_feature_set(d12)
x413 = generate_feature_set(d13)
x414 = generate_feature_set(d14)
x415 = generate_feature_set(d15)
x416 = generate_feature_set(d16)

print " Roll read done"


x51 = generate_feature_set_kmer(e1)
x52 = generate_feature_set_kmer(e2)
x53 = generate_feature_set_kmer(e3)
x54 = generate_feature_set_kmer(e4)
x55 = generate_feature_set_kmer(e5)
x56 = generate_feature_set_kmer(e6)
x57 = generate_feature_set_kmer(e7)
x58 = generate_feature_set_kmer(e8)
x59 = generate_feature_set_kmer(e9)
x510 = generate_feature_set_kmer(e10)
x511 = generate_feature_set_kmer(e11)
x512 = generate_feature_set_kmer(e12)
x513 = generate_feature_set_kmer(e13)
x514 = generate_feature_set_kmer(e14)
x515 = generate_feature_set_kmer(e15)
x516 = generate_feature_set_kmer(e16)

print "kmer read done"


M = np.zeros([4,20])
M[0,0] = 0.247465
M[0,1] = 0.237323
M[0,2] = 0.141988
M[0,3] = 0.093306
M[0,4] = 0.131846
M[0,5] = 0.338742
M[0,6] = 0.442191
M[0,7] = 0.403651
M[0,8] = 0.141988
M[0,9] = 0.046653
M[0,10] = 0.028398
M[0,11] = 0.010142
M[0,12] = 0.000000
M[0,13] = 0.002028
M[0,14] = 0.000000
M[0,15] = 0.000000
M[0,16] = 0.008114
M[0,17] = 0.265720
M[0,18] = 0.229209
M[0,19] = 0.135903

M[1,0] = 0.338742
M[1,1] = 0.257606
M[1,2] = 0.269777
M[1,3] = 0.298174
M[1,4] = 0.387424
M[1,5] = 0.283976
M[1,6] = 0.217039
M[1,7] = 0.170385
M[1,8] = 0.273834
M[1,9] = 0.263692
M[1,10] = 0.103448
M[1,11] = 0.795132
M[1,12] = 0.679513
M[1,13] = 0.995943
M[1,14] = 0.000000
M[1,15] = 0.937120
M[1,16] = 0.610548
M[1,17] = 0.359026
M[1,18] = 0.306288
M[1,19] = 0.430020

M[2,0] = 0.233266
M[2,1] = 0.237323
M[2,2] = 0.180527
M[2,3] = 0.174442
M[2,4] = 0.198783
M[2,5] = 0.298174
M[2,6] = 0.186613
M[2,7] = 0.158215
M[2,8] = 0.105477
M[2,9] = 0.034483
M[2,10] = 0.115619
M[2,11] = 0.176471
M[2,12] = 0.320487
M[2,13] = 0.000000
M[2,14] = 0.983773
M[2,15] = 0.062880
M[2,16] = 0.286004
M[2,17] = 0.180527
M[2,18] = 0.257606
M[2,19] = 0.308316

M[3,0] = 0.180527
M[3,1] = 0.267748
M[3,2] = 0.407708
M[3,3] = 0.434077
M[3,4] = 0.281947
M[3,5] = 0.079108
M[3,6] = 0.154158
M[3,7] = 0.267748
M[3,8] = 0.478702
M[3,9] = 0.655172
M[3,10] = 0.752535
M[3,11] = 0.018256
M[3,12] = 0.000000
M[3,13] = 0.002028
M[3,14] = 0.016227
M[3,15] = 0.000000
M[3,16] = 0.095335
M[3,17] = 0.194726
M[3,18] = 0.206897
M[3,19] = 0.125761


x61 = get_pwm_score(e1,M)
x62 = get_pwm_score(e2,M)
x63 = get_pwm_score(e3,M)
x64 = get_pwm_score(e4,M)
x65 = get_pwm_score(e5,M)
x66 = get_pwm_score(e6,M)
x67 = get_pwm_score(e7,M)
x68 = get_pwm_score(e8,M)
x69 = get_pwm_score(e9,M)
x610 = get_pwm_score(e10,M)
x611 = get_pwm_score(e11,M)
x612 = get_pwm_score(e12,M)
x613 = get_pwm_score(e13,M)
x614 = get_pwm_score(e14,M)
x615 = get_pwm_score(e15,M)
x616 = get_pwm_score(e16,M)

print "PWM read done"

Xpos1 = np.vstack([x11,x12,x13,x14,x15,x16,x17,x18,x19,x110,x111,x112,x113,x114,x115,x116])
Xpos2 = np.vstack([x21,x22,x23,x24,x25,x26,x27,x28,x29,x210,x211,x212,x213,x214,x215,x216])
Xpos3 = np.vstack([x31,x32,x33,x34,x35,x36,x37,x38,x39,x310,x311,x312,x313,x314,x315,x316])
Xpos4 = np.vstack([x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416])
Xpos5 = np.vstack([x51,x52,x53,x54,x55,x56,x57,x58,x59,x510,x511,x512,x513,x514,x515,x516])
Xpos6 = np.vstack([x61,x62,x63,x64,x65,x66,x67,x68,x69,x610,x611,x612,x613,x614,x615,x616])
Xpos = np.hstack([Xpos1,Xpos2,Xpos3,Xpos4,Xpos5,Xpos6])
pickle.dump(Xpos,open('positiveFeatureSet_1-22_E2F4.p',"wb"))






