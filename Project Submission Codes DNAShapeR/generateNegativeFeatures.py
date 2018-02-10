import itertools 
from sklearn import preprocessing
import numpy as np
import random
import pickle
import re
import math
import pickle

def generate_feature_set(shapelist):
	x = np.zeros(20)
	for k,j in enumerate(shapelist):
		if(k%1000==0 and k!=0):
			print k,'of',len(shapelist)
		x = np.vstack([x,j])
	x = x[1:]
	return x


def get_pwm_score(shapelist,M):
	x = np.zeros(1)
	for l,i in enumerate(shapelist):
		#print i
		p = 1
		if (l%1000==0 and l!=0):
			print l,'of',len(shapelist)
			#break
		for j,k in enumerate(i):
			#print k,p
			if(k=='A'):
				print M[0,j]
				p*=M[0,j]
			elif(k=='C'):
				print M[1,j]
				p*=M[1,j]
			elif(k=='G'):
				print M[2,j]
				p*=M[2,j]
			elif(k=='T'):
				print M[3,j]
				p*=M[3,j]
			else:
				p*=0
		x = np.vstack([x,p])
	x = x[1:]
	return x

def generate_feature_set_kmer(shapelist):
	x = np.zeros(20)
	for z,i in enumerate(shapelist):
		if (z%1000==0 and z!=0):
			print z,'of',len(shapelist)
			#break
		#print i
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
				
		#print p
		#print x
		
		x = np.vstack([x,p])
	x = x[1:]
	return x

a29 = pickle.load(open('featuresetNegative/negativeShapeFeatureHeltchr20.p',"rb"))

a30 = pickle.load(open('featuresetNegative/negativeShapeFeatureHeltchr22.p',"rb"))
print 'loaded helt pickle files'


b29 = pickle.load(open('featuresetNegative/negativeShapeFeatureMgwchr20.p',"rb"))

b30 = pickle.load(open('featuresetNegative/negativeShapeFeatureMgwchr22.p',"rb"))

print 'loaded mgw pickle files'

c29 = pickle.load(open('featuresetNegative/negativeShapeFeatureProtchr20.p',"rb"))

c30 = pickle.load(open('featuresetNegative/negativeShapeFeatureProtchr22.p',"rb"))
print 'loaded prot pickle files'

d29 = pickle.load(open('featuresetNegative/negativeShapeFeatureRollchr20.p',"rb"))

d30 = pickle.load(open('featuresetNegative/negativeShapeFeatureRollchr22.p',"rb"))

print 'loaded roll pickle files'

e29 = pickle.load(open('featuresetNegative/negativeKmerFeaturechr20.p',"rb"))

e30 = pickle.load(open('featuresetNegative/negativeKmerFeaturechr22.p',"rb"))

print 'loaded kmer pickle files'
print 'Started helt'
x129 = generate_feature_set(a29)
print 'Done 13'

x130 = generate_feature_set(a30)
print 'Done 14'
print " Helt read done"

print 'started reading mgw'
x229 = generate_feature_set(b29)
print 'Done 13'

x230 = generate_feature_set(b30)
print 'Done 14'
print " MGW read done"

print 'started prot'
x329 = generate_feature_set(c29)
print 'Done 13'

x330 = generate_feature_set(c30)
print 'Done 14'
print " prot read done"

print 'started roll'
x429 = generate_feature_set(d29)
print 'Done 13'

x430 = generate_feature_set(d30)
print 'Done 14'

print " Roll read done"
print 'started kmer'
x529 = generate_feature_set_kmer(e29)
print e29[0]
print 'Done 13'

x530 = generate_feature_set_kmer(e30)
print 'Done 14'
print "kmer read done"
'''
#=====================================Pwm matrix for Max========================================
M = np.zeros([4,14])
M[0,0] = 0.210953
M[0,1] = 0.391481
M[0,2] = 0.146045
M[0,3] = 0.000000
M[0,4] = 0.997972
M[0,5] = 0.000000
M[0,6] = 0.000000
M[0,7] = 0.000000
M[0,8] = 0.000000
M[0,9] = 0.336714
M[0,10] = 0.020284
M[0,11] = 0.127789
M[0,12] = 0.186613
M[0,13] = 0.101420
#0.263692,0.247465,0.588235,1.000000,0.000000,0.929006,0.000000,0.004057,0.000000,0.231237,0.527383,0.464503,0.320487,0.381339,
M[1,0] = 0.263692
M[1,1] = 0.247465
M[1,2] = 0.588235
M[1,3] = 1.000000
M[1,4] = 0.000000
M[1,5] = 0.929006
M[1,6] = 0.000000
M[1,7] = 0.004057
M[1,8] = 0.000000
M[1,9] = 0.231237
M[1,10] = 0.527383
M[1,11] = 0.464503
M[1,12] = 0.320487
M[1,13] = 0.381339
#0.407708,0.219067,0.229209,0.000000,0.002028,0.000000,1.000000,0.000000,1.000000,0.249493,0.101420,0.052738,0.208925,0.310345,
M[2,0] = 0.407708
M[2,1] = 0.219067
M[2,2] = 0.229209
M[2,3] = 0.000000
M[2,4] = 0.002028
M[2,5] = 0.000000
M[2,6] = 1.000000
M[2,7] = 0.000000
M[2,8] = 1.000000
M[2,9] = 0.249493
M[2,10] = 0.101420
M[2,11] = 0.052738
M[2,12] = 0.208925
M[2,13] = 0.310345
0.117647,0.141988,0.036511,0.000000,0.000000,0.070994,0.000000,0.995943,0.000000,0.182556,0.350913,0.354970,0.283976,0.206897,
M[3,0] = 0.117647
M[3,1] = 0.141988
M[3,2] = 0.036511
M[3,3] = 0.000000
M[3,4] = 0.000000
M[3,5] = 0.070994
M[3,6] = 0.000000
M[3,7] = 0.995943
M[3,8] = 0.000000
M[3,9] = 0.182556
M[3,10] = 0.350913
M[3,11] = 0.354970
M[3,12] = 0.283976
M[3,13] = 0.206897
'''
#=====================================Pwm matrix for ZNF263=====================================
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
print 'started pwm'
x629 = get_pwm_score(e29,M)

x630 = get_pwm_score(e30,M)

print "PWM read done"

Xneg1 = np.vstack([x129,x130])

Xneg2 = np.vstack([x229,x230])
Xneg3 = np.vstack([x329,x330])
Xneg4 = np.vstack([x429,x430])
Xneg5 = np.vstack([x529,x530])
Xneg6 = np.vstack([x629,x630])
Xneg = np.hstack([Xneg1,Xneg2,Xneg3,Xneg4,Xneg5,Xneg6])
pickle.dump(Xneg,open('negativeFeatureSet_21_22.p',"wb"))
