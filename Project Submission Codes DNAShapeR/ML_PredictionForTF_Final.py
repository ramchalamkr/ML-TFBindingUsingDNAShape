# ML algorithms - prediction of binding sites. alogn with the classifier scores.
import itertools 
from sklearn import preprocessing
import matplotlib.pyplot as plt
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
from sklearn import svm
from sklearn.metrics import precision_recall_curve
from itertools import cycle
from sklearn.metrics import confusion_matrix
from sklearn.metrics import roc_curve

def TN_TP(ytrue,ypred,l,tpe):
	count = 0
	for i in range(l):
		if(ytrue[i] == tpe and ypred[i] == tpe):
			count+=1
	return count

def FN_FP(ytrue,ypred,l,tpe):
	count = 0
	for i in range(l):
		if(ytrue[i] == tpe and ypred[i] != tpe):
			count+=1
	return count

def calculate_TPR_FPR(ytrue,ypred,P):
	#T1 = acuracy(ytrue,ypred,len(P))
	#print "accuracy:" + str(T1)

	TP = TN_TP(ytrue,ypred,len(P),1)
	TN = TN_TP(ytrue,ypred,len(P),0)
	FP = FN_FP(ytrue,ypred,len(P),0)
	FN = FN_FP(ytrue,ypred,len(P),1)
	print "False positivies"
	print FP
	print "False Negatives"
	print FN
	TPR = float(TP)/(TP+FN)
	FPR = float(FP)/(FP+TN)

	return TPR,FPR


Xneg = pickle.load(open('negativeFeatureSet_20_22.p',"rb"))
#Xpos = pickle.load(open('positiveFeatureSet_1-22_E2F4.p',"rb"))
#Xpos = pickle.load(open('positiveFeatureSet_1-22_MAX.p',"rb"))
Xpos = pickle.load(open('positiveFeatureSet_1-22_ZNF263.p',"rb"))
#Xpos = pickle.load(open('positiveFeatureSet_1-22_UAK21.p',"rb"))

'''
#============================For 14 length MAX=========================
XnegHelt = Xneg[:,0:14]
#print XnegHelt[0]
#print len(XnegHelt)
XnegMgw = Xneg[:,20:34]
#print XnegMgw[0]
#print len(XnegMgw)
XnegProt = Xneg[:,40:54]
#print XnegProt[0]
#print len(XnegProt)
XnegRoll = Xneg[:,60:74]
#print XnegRoll[0]
#print len(XnegRoll)
XnegKmer = Xneg[:,80:94]
#print XnegKmer[0]
#print len(XnegKmer)
Xpwm = Xneg[:,100]
z=[]
for i in Xpwm:
	z.append([i])
#print len(z)
Xneg = np.hstack([XnegHelt,XnegMgw,XnegProt,XnegRoll,XnegKmer,z])
#==================================================================
'''


#============================For 15 length ZNF263=========================
XnegHelt = Xneg[:,0:15]
#print XnegHelt[0]
#print len(XnegHelt)
XnegMgw = Xneg[:,20:35]
#print XnegMgw[0]
#print len(XnegMgw)
XnegProt = Xneg[:,40:55]
#print XnegProt[0]
#print len(XnegProt)
XnegRoll = Xneg[:,60:75]
#print XnegRoll[0]
#print len(XnegRoll)
XnegKmer = Xneg[:,80:95]
#print XnegKmer[0]
#print len(XnegKmer)
Xpwm = Xneg[:,100]
z=[]
for i in Xpwm:
	z.append([i])
#print len(z)
Xneg = np.hstack([XnegHelt,XnegMgw,XnegProt,XnegRoll,XnegKmer,z])
#==================================================================

Xneg = Xneg[:,60:75]
Xpos = Xpos[:,60:75]

Xpos = [s for s in Xpos if "NA" not in s and "abcd" not in s]
Ypos = [1 for x in range(len(Xpos))]

Xneg = [s for s in Xneg if "NA" not in s and "abcd" not in s]
Yneg = [0 for x in range(len(Xneg))]

X0 = np.vstack([Xpos,Xneg])
Y = Ypos+Yneg
X = X0.astype(np.float)

print "no of positive samples"
print len(Ypos)
print "No of negative samples"
print len(Yneg)

X_train, X_test, y_train, y_test = train_test_split(X,Y,test_size=0.33, random_state=42)
print "Training Set count"
print len(X_train)
print "Test Set Count"
print len(X_test)
#clf = svm.SVC()
#clf = RandomForestClassifier(n_estimators=20)
clf = ExtraTreesClassifier(n_estimators=20)
#clf = MultinomialNB()
#clf = MLPClassifier(solver ='adam',alpha=1e-6,hidden_layer_sizes = (500,),random_state = 1)
clf = clf.fit(X_train, y_train)
y_pred = clf.predict(X_test)
print "accuracy"
print np.mean(y_pred == y_test)
print "error count"
print np.sum(y_pred != y_test)
c=0
for i in range(len(y_test)):
	if(y_test[i] == 0 and y_pred[i] == 1):
		c+=1

print "false positive" + str(c)
print "no of examples negative"
print len([x for x in y_pred if x ==0])
print "no of examples positive"
print len([x for x in y_pred if x ==1])

precision, recall, thresholds = precision_recall_curve(y_test,y_pred)
'''
colors = cycle(['navy', 'turquoise', 'darkorange', 'cornflowerblue', 'teal'])
lw = 2
plt.clf()
plt.plot(recall, precision, lw=lw, color='navy',
         label='Precision-Recall curve')
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.05])
plt.xlim([0.0, 1.05])
plt.title('Precision-Recall example: UAK21 - Feature Set D')
plt.legend(loc="lower left")
plt.show()
'''

TPR,FPR = calculate_TPR_FPR(y_test,y_pred,y_test)
print "sensitivity"
print TPR
print "specificity"
print 1-FPR
print "10 fold CV accuracy"
#from sklearn.metrics import precision_recall_curve

#predicted = cross_val_predict(clf,X ,Y, cv=10)
#print metrics.accuracy_score(Y, predicted)
#target_names = ['class NTF', 'class TF']
#print(classification_report(y_test, y_pred, target_names=target_names))
