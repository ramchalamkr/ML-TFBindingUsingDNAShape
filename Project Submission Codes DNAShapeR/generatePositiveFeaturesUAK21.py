#To generate the positive features for UAK21 TF
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
	for k,j in enumerate(shapelist):
		#print j
		x = np.vstack([x,j])
		if(k%100==0 and k!=0):
			#break
			print k,'of',len(shapelist)
		
	x = x[1:]
	return x


def get_pwm_score(shapelist,M):
	x = np.zeros(1)
	for l,i in enumerate(shapelist):
		p = 1
		if (l%100==0 and l!=0):
			print l,'of',len(shapelist)
			#break
		for j,k in enumerate(i):
			#print k,p
			#if(j==15):
			#	continue
			if(k=='A' or k=='a'):
				#print M[0,j]
				p*=M[0,j]
			elif(k=='C' or k=='c'):
				#print M[1,j]
				p*=M[1,j]
			elif(k=='G' or k=='g'):
				#print M[2,j]
				p*=M[2,j]
			elif(k=='T' or k=='t'):
				#print M[3,j]
				p*=M[3,j]
			else:
				p*=0
		x = np.vstack([x,p])
	x = x[1:]
	return x

def generate_feature_set_kmer(shapelist):
	x = np.zeros(20)
	for z,i in enumerate(shapelist):
		if (z%100==0 and z!=0):
			print z,'of',len(shapelist)
		p = []
		for a,k in enumerate(i):
			#print k
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

a1 = pickle.load(open('UAK21/UAK21_helt_Positivechr1.p',"rb"))
a2 = pickle.load(open('UAK21/UAK21_helt_Positivechr2.p',"rb"))
a3 = pickle.load(open('UAK21/UAK21_helt_Positivechr3.p',"rb"))
a4 = pickle.load(open('UAK21/UAK21_helt_Positivechr4.p',"rb"))
a5 = pickle.load(open('UAK21/UAK21_helt_Positivechr5.p',"rb"))
a6 = pickle.load(open('UAK21/UAK21_helt_Positivechr6.p',"rb"))
a7 = pickle.load(open('UAK21/UAK21_helt_Positivechr7.p',"rb"))
a8 = pickle.load(open('UAK21/UAK21_helt_Positivechr8.p',"rb"))
a9 = pickle.load(open('UAK21/UAK21_helt_Positivechr9.p',"rb"))
a10 = pickle.load(open('UAK21/UAK21_helt_Positivechr10.p',"rb"))
a11 = pickle.load(open('UAK21/UAK21_helt_Positivechr11.p',"rb"))
a12 = pickle.load(open('UAK21/UAK21_helt_Positivechr12.p',"rb"))
a13 = pickle.load(open('UAK21/UAK21_helt_Positivechr13.p',"rb"))
a14 = pickle.load(open('UAK21/UAK21_helt_Positivechr14.p',"rb"))
a15 = pickle.load(open('UAK21/UAK21_helt_Positivechr15.p',"rb"))
a16 = pickle.load(open('UAK21/UAK21_helt_Positivechr16.p',"rb"))
a17 = pickle.load(open('UAK21/UAK21_helt_Positivechr17.p',"rb"))
a18 = pickle.load(open('UAK21/UAK21_helt_Positivechr18.p',"rb"))
a19 = pickle.load(open('UAK21/UAK21_helt_Positivechr19.p',"rb"))
a20 = pickle.load(open('UAK21/UAK21_helt_Positivechr20.p',"rb"))
a21 = pickle.load(open('UAK21/UAK21_helt_Positivechr21.p',"rb"))
a22 = pickle.load(open('UAK21/UAK21_helt_Positivechr22.p',"rb"))
a23 = pickle.load(open('UAK21/UAK21_helt_PositivechrX.p',"rb"))

print 'loaded helt pickle files'
b1 = pickle.load(open('UAK21/UAK21_mgw_Positivechr1.p',"rb"))
b2 = pickle.load(open('UAK21/UAK21_mgw_Positivechr2.p',"rb"))
b3 = pickle.load(open('UAK21/UAK21_mgw_Positivechr3.p',"rb"))
b4 = pickle.load(open('UAK21/UAK21_mgw_Positivechr4.p',"rb"))
b5 = pickle.load(open('UAK21/UAK21_mgw_Positivechr5.p',"rb"))
b6 = pickle.load(open('UAK21/UAK21_mgw_Positivechr6.p',"rb"))
b7 = pickle.load(open('UAK21/UAK21_mgw_Positivechr7.p',"rb"))
b8 = pickle.load(open('UAK21/UAK21_mgw_Positivechr8.p',"rb"))
b9 = pickle.load(open('UAK21/UAK21_mgw_Positivechr9.p',"rb"))
b10 = pickle.load(open('UAK21/UAK21_mgw_Positivechr10.p',"rb"))
b11 = pickle.load(open('UAK21/UAK21_mgw_Positivechr11.p',"rb"))
b12 = pickle.load(open('UAK21/UAK21_mgw_Positivechr12.p',"rb"))
b13 = pickle.load(open('UAK21/UAK21_mgw_Positivechr13.p',"rb"))
b14 = pickle.load(open('UAK21/UAK21_mgw_Positivechr14.p',"rb"))
b15 = pickle.load(open('UAK21/UAK21_mgw_Positivechr15.p',"rb"))
b16 = pickle.load(open('UAK21/UAK21_mgw_Positivechr16.p',"rb"))
b17 = pickle.load(open('UAK21/UAK21_mgw_Positivechr17.p',"rb"))
b18 = pickle.load(open('UAK21/UAK21_mgw_Positivechr18.p',"rb"))
b19 = pickle.load(open('UAK21/UAK21_mgw_Positivechr19.p',"rb"))
b20 = pickle.load(open('UAK21/UAK21_mgw_Positivechr20.p',"rb"))
b21 = pickle.load(open('UAK21/UAK21_mgw_Positivechr21.p',"rb"))
b22 = pickle.load(open('UAK21/UAK21_mgw_Positivechr22.p',"rb"))
b23 = pickle.load(open('UAK21/UAK21_mgw_PositivechrX.p',"rb"))


print 'loaded mgw pickle files'

c1 = pickle.load(open('UAK21/UAK21_prot_Positivechr1.p',"rb"))
c2 = pickle.load(open('UAK21/UAK21_prot_Positivechr2.p',"rb"))
c3 = pickle.load(open('UAK21/UAK21_prot_Positivechr3.p',"rb"))
c4 = pickle.load(open('UAK21/UAK21_prot_Positivechr4.p',"rb"))
c5 = pickle.load(open('UAK21/UAK21_prot_Positivechr5.p',"rb"))
c6 = pickle.load(open('UAK21/UAK21_prot_Positivechr6.p',"rb"))
c7 = pickle.load(open('UAK21/UAK21_prot_Positivechr7.p',"rb"))
c8 = pickle.load(open('UAK21/UAK21_prot_Positivechr8.p',"rb"))
c9 = pickle.load(open('UAK21/UAK21_prot_Positivechr9.p',"rb"))
c10 = pickle.load(open('UAK21/UAK21_prot_Positivechr10.p',"rb"))
c11 = pickle.load(open('UAK21/UAK21_prot_Positivechr11.p',"rb"))
c12 = pickle.load(open('UAK21/UAK21_prot_Positivechr12.p',"rb"))
c13 = pickle.load(open('UAK21/UAK21_prot_Positivechr13.p',"rb"))
c14 = pickle.load(open('UAK21/UAK21_prot_Positivechr14.p',"rb"))
c15 = pickle.load(open('UAK21/UAK21_prot_Positivechr15.p',"rb"))
c16 = pickle.load(open('UAK21/UAK21_prot_Positivechr16.p',"rb"))
c17 = pickle.load(open('UAK21/UAK21_prot_Positivechr17.p',"rb"))
c18 = pickle.load(open('UAK21/UAK21_prot_Positivechr18.p',"rb"))
c19 = pickle.load(open('UAK21/UAK21_prot_Positivechr19.p',"rb"))
c20 = pickle.load(open('UAK21/UAK21_prot_Positivechr20.p',"rb"))
c21 = pickle.load(open('UAK21/UAK21_prot_Positivechr21.p',"rb"))
c22 = pickle.load(open('UAK21/UAK21_prot_Positivechr22.p',"rb"))
c23 = pickle.load(open('UAK21/UAK21_prot_PositivechrX.p',"rb"))

print 'loaded prot pickle files'

d1 = pickle.load(open('UAK21/UAK21_roll_Positivechr1.p',"rb"))
d2 = pickle.load(open('UAK21/UAK21_roll_Positivechr2.p',"rb"))
d3 = pickle.load(open('UAK21/UAK21_roll_Positivechr3.p',"rb"))
d4 = pickle.load(open('UAK21/UAK21_roll_Positivechr4.p',"rb"))
d5 = pickle.load(open('UAK21/UAK21_roll_Positivechr5.p',"rb"))
d6 = pickle.load(open('UAK21/UAK21_roll_Positivechr6.p',"rb"))
d7 = pickle.load(open('UAK21/UAK21_roll_Positivechr7.p',"rb"))
d8 = pickle.load(open('UAK21/UAK21_roll_Positivechr8.p',"rb"))
d9 = pickle.load(open('UAK21/UAK21_roll_Positivechr9.p',"rb"))
d10 = pickle.load(open('UAK21/UAK21_roll_Positivechr10.p',"rb"))
d11 = pickle.load(open('UAK21/UAK21_roll_Positivechr11.p',"rb"))
d12 = pickle.load(open('UAK21/UAK21_roll_Positivechr12.p',"rb"))
d13 = pickle.load(open('UAK21/UAK21_roll_Positivechr13.p',"rb"))
d14 = pickle.load(open('UAK21/UAK21_roll_Positivechr14.p',"rb"))
d15 = pickle.load(open('UAK21/UAK21_roll_Positivechr15.p',"rb"))
d16 = pickle.load(open('UAK21/UAK21_roll_Positivechr16.p',"rb"))
d17 = pickle.load(open('UAK21/UAK21_roll_Positivechr17.p',"rb"))
d18 = pickle.load(open('UAK21/UAK21_roll_Positivechr18.p',"rb"))
d19 = pickle.load(open('UAK21/UAK21_roll_Positivechr19.p',"rb"))
d20 = pickle.load(open('UAK21/UAK21_roll_Positivechr20.p',"rb"))
d21 = pickle.load(open('UAK21/UAK21_roll_Positivechr21.p',"rb"))
d22 = pickle.load(open('UAK21/UAK21_roll_Positivechr22.p',"rb"))
d23 = pickle.load(open('UAK21/UAK21_roll_PositivechrX.p',"rb"))


print 'loaded roll pickle files'

e1 = pickle.load(open('UAK21/UAK21_kmer_Positivechr1.p',"rb"))
e2 = pickle.load(open('UAK21/UAK21_kmer_Positivechr2.p',"rb"))
e3 = pickle.load(open('UAK21/UAK21_kmer_Positivechr3.p',"rb"))
e4 = pickle.load(open('UAK21/UAK21_kmer_Positivechr4.p',"rb"))
e5 = pickle.load(open('UAK21/UAK21_kmer_Positivechr5.p',"rb"))
e6 = pickle.load(open('UAK21/UAK21_kmer_Positivechr6.p',"rb"))
e7 = pickle.load(open('UAK21/UAK21_kmer_Positivechr7.p',"rb"))
e8 = pickle.load(open('UAK21/UAK21_kmer_Positivechr8.p',"rb"))
e9 = pickle.load(open('UAK21/UAK21_kmer_Positivechr9.p',"rb"))
e10 = pickle.load(open('UAK21/UAK21_kmer_Positivechr10.p',"rb"))
e11 = pickle.load(open('UAK21/UAK21_kmer_Positivechr11.p',"rb"))
e12 = pickle.load(open('UAK21/UAK21_kmer_Positivechr12.p',"rb"))
e13 = pickle.load(open('UAK21/UAK21_kmer_Positivechr13.p',"rb"))
e14 = pickle.load(open('UAK21/UAK21_kmer_Positivechr14.p',"rb"))
e15 = pickle.load(open('UAK21/UAK21_kmer_Positivechr15.p',"rb"))
e16 = pickle.load(open('UAK21/UAK21_kmer_Positivechr16.p',"rb"))
e17 = pickle.load(open('UAK21/UAK21_kmer_Positivechr17.p',"rb"))
e18 = pickle.load(open('UAK21/UAK21_kmer_Positivechr18.p',"rb"))
e19 = pickle.load(open('UAK21/UAK21_kmer_Positivechr19.p',"rb"))
e20 = pickle.load(open('UAK21/UAK21_kmer_Positivechr20.p',"rb"))
e21 = pickle.load(open('UAK21/UAK21_kmer_Positivechr21.p',"rb"))
e22 = pickle.load(open('UAK21/UAK21_kmer_Positivechr22.p',"rb"))
e23 = pickle.load(open('UAK21/UAK21_kmer_PositivechrX.p',"rb"))

print 'loaded kmer pickle files'

x11 = generate_feature_set(a1)
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
x117 = generate_feature_set(a17)
x118 = generate_feature_set(a18)
x119 = generate_feature_set(a19)
x120 = generate_feature_set(a20)
x121 = generate_feature_set(a21)
x122 = generate_feature_set(a22)
x123 = generate_feature_set(a23)


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
x217 = generate_feature_set(b17)
x218 = generate_feature_set(b18)
x219 = generate_feature_set(b19)
x220 = generate_feature_set(b20)
x221 = generate_feature_set(b21)
x222 = generate_feature_set(b22)
x223 = generate_feature_set(b23)

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
x317 = generate_feature_set(c17)
x318 = generate_feature_set(c18)
x319 = generate_feature_set(c19)
x320 = generate_feature_set(c20)
x321 = generate_feature_set(c21)
x322 = generate_feature_set(c22)
x323 = generate_feature_set(c23)


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
x417 = generate_feature_set(d17)
x418 = generate_feature_set(d18)
x419 = generate_feature_set(d19)
x420 = generate_feature_set(d20)
x421 = generate_feature_set(d21)
x422 = generate_feature_set(d22)
x423 = generate_feature_set(d23)

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
x517 = generate_feature_set_kmer(e17)
x518 = generate_feature_set_kmer(e18)
x519 = generate_feature_set_kmer(e19)
x520 = generate_feature_set_kmer(e20)
x521 = generate_feature_set_kmer(e21)
x522 = generate_feature_set_kmer(e22)
x523 = generate_feature_set_kmer(e23)

print "kmer read done"


M = np.zeros([4,20])
#0.000000,0.000000,0.000000,1.000000,0.971014,0.028986,0.000000,0.144928,0.275362,0.275362,
#0.115942,0.144928,0.144928,0.072464,0.101449,0.159420,0.144928,0.086957,0.144928,0.101449,
M[0,0] = 0.000000
M[0,1] = 0.000000
M[0,2] = 0.000000
M[0,3] = 1.000000
M[0,4] = 0.971014
M[0,5] = 0.028986
M[0,6] = 0.000000
M[0,7] = 0.144928
M[0,8] = 0.275362
M[0,9] = 0.275362
M[0,10] = 0.115942
M[0,11] = 0.144928
M[0,12] = 0.144928
M[0,13] = 0.072464
M[0,14] = 0.101449
M[0,15] = 0.159420
M[0,16] = 0.144928
M[0,17] = 0.086957
M[0,18] = 0.144928
M[0,19] = 0.101449
#0.927536,0.000000,0.000000,0.000000,0.000000,0.971014,0.739130,0.318841,0.072464,0.115942,
#0.275362,0.159420,0.217391,0.608696,0.275362,0.376812,0.362319,0.159420,0.115942,0.536232,
M[1,0] = 0.927536
M[1,1] = 0.000000
M[1,2] = 0.000000
M[1,3] = 0.000000
M[1,4] = 0.000000
M[1,5] = 0.971014
M[1,6] = 0.739130
M[1,7] = 0.318841
M[1,8] = 0.072464
M[1,9] = 0.115942
M[1,10] = 0.275362
M[1,11] = 0.159420
M[1,12] = 0.217391
M[1,13] = 0.608696
M[1,14] = 0.275362
M[1,15] = 0.376812
M[1,16] = 0.362319
M[1,17] = 0.159420
M[1,18] = 0.115942
M[1,19] = 0.536232
#0.072464,1.000000,1.000000,0.000000,0.028986,0.000000,0.130435,0.028986,0.217391,0.260870,
#0.550725,0.260870,0.434783,0.130435,0.347826,0.391304,0.318841,0.739130,0.608696,0.289855,
M[2,0] = 0.072464
M[2,1] = 1.000000
M[2,2] = 1.000000
M[2,3] = 0.000000
M[2,4] = 0.028986
M[2,5] = 0.000000
M[2,6] = 0.130435
M[2,7] = 0.028986
M[2,8] = 0.217391
M[2,9] = 0.260870
M[2,10] = 0.550725
M[2,11] = 0.260870
M[2,12] = 0.434783
M[2,13] = 0.130435
M[2,14] = 0.347826
M[2,15] = 0.391304
M[2,16] = 0.318841
M[2,17] = 0.739130
M[2,18] = 0.608696
M[2,19] = 0.289855
#0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.130435,0.507246,0.434783,0.347826,
#0.057971,0.434783,0.202899,0.188406,0.275362,0.072464,0.173913,0.014493,0.130435,0.072464,
M[3,0] = 0.000000
M[3,1] = 0.000000
M[3,2] = 0.000000
M[3,3] = 0.000000
M[3,4] = 0.000000
M[3,5] = 0.000000
M[3,6] = 0.130435
M[3,7] = 0.507246
M[3,8] = 0.434783
M[3,9] = 0.347826
M[3,10] = 0.057971
M[3,11] = 0.434783
M[3,12] = 0.202899
M[3,13] = 0.188406
M[3,14] = 0.275362
M[3,15] = 0.072464
M[3,16] = 0.173913
M[3,17] = 0.014493
M[3,18] = 0.130435
M[3,19] = 0.072464



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
x617 = get_pwm_score(e17,M)
x618 = get_pwm_score(e18,M)
x619 = get_pwm_score(e19,M)
x620 = get_pwm_score(e20,M)
x621 = get_pwm_score(e21,M)
x622 = get_pwm_score(e22,M)
x623 = get_pwm_score(e23,M)

print "PWM read done"

Xpos1 = np.vstack([x11,x12,x13,x14,x15,x16,x17,x18,x19,x110,x111,x112,x113,x114,x115,x116,x117,x118,x119,x120,x121,x122,x123])
Xpos2 = np.vstack([x21,x22,x23,x24,x25,x26,x27,x28,x29,x210,x211,x212,x213,x214,x215,x216,x217,x218,x219,x220,x221,x222,x223])
Xpos3 = np.vstack([x31,x32,x33,x34,x35,x36,x37,x38,x39,x310,x311,x312,x313,x314,x315,x316,x317,x318,x319,x320,x321,x322,x323])
Xpos4 = np.vstack([x41,x42,x43,x44,x45,x46,x47,x48,x49,x410,x411,x412,x413,x414,x415,x416,x417,x418,x419,x420,x421,x422,x423])
Xpos5 = np.vstack([x51,x52,x53,x54,x55,x56,x57,x58,x59,x510,x511,x512,x513,x514,x515,x516,x517,x518,x519,x520,x521,x522,x523])
Xpos6 = np.vstack([x61,x62,x63,x64,x65,x66,x67,x68,x69,x610,x611,x612,x613,x614,x615,x616,x617,x618,x619,x620,x621,x622,x623])
Xpos = np.hstack([Xpos1,Xpos2,Xpos3,Xpos4,Xpos5,Xpos6])
pickle.dump(Xpos,open('positiveFeatureSet_1-22_UAK21.p',"wb"))




