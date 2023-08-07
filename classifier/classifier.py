from __future__ import print_function
import datetime
from isoelectric import ipc #Requires Isoelectric Package
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import statistics
from propy import PyPro #Requires Installing Propy Package
import propy #Requires Installing Propy Package
from propy.PyPro import GetProDes #Requires Installing Propy Package
import peptides
import hmomentk

import pickle


import argparse
import math
import os
import time
labelDict={'negative':'non-AMP' , 'positive': 'AMP'}

aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}

#
b= open('input_peptide_sequences.txt').read().splitlines()
def Merge(dict1, dict2): 
    return(dict2.update(dict1))
pp=[]

#Create Dictionary of Sequences
c={}
x = range(1,len(b)+1)
dic={}
for i in b:
    dic[b.index(i)+1]=i
#Calculate Properties
aggDict = pickle.load(open('aggDict.pkl', 'rb'))#Requires aggDict File
chargeDict = pickle.load(open('chargeDict.pkl', 'rb'))#Requires chargeDict File
#isoElP=[]
netcharge=[]
agg=[]
chargeDensities=[]
mw=[]
hm=[]
hp=[]
isop=[]

dsss={"seq":b}
for key, val in dic.items():
    agList=[]
    chargeList=[]
    mwList=[]
    for j in list(val):
        agList.append(aggDict[j])
        chargeList.append(chargeDict[j])
    agg.append(statistics.mean(agList))
    chargeDensities.append(sum(chargeList)/ipc.calculate_molecular_weight(i))
    mw.append(ipc.calculate_molecular_weight(val))
    netcharge.append(sum(chargeList))
    pepseq = peptides.Peptide(val)
    hp.append(-1*peptides.Peptide.hydrophobicity(pepseq,scale="KyteDoolittle"))
    hm.append(hmomentk.hydrop_moment(pepseq))
    isop.append(peptides.Peptide.isoelectric_point(pepseq,pKscale="Lehninger"))
    
dsss["molecularWeight"]=mw
dsss["netCharge"]=netcharge
dsss["chargeDensity"]=chargeDensities
dsss["aggregationPropensityInVivo"]=agg
dsss["hydrophobicity"]=hp
dsss["hydrophobicMoment"]=hm
dsss["isoelectricPoint"]=isop
df=pd.DataFrame(data=pp)
p5=pd.DataFrame(data=dsss)
p6=pd.concat([p5, df.reindex(p5.index)], axis=1)
#Rescaling Data
#Needs Final Features MinMaxes
p=[]
features=["molecularWeight","netCharge","chargeDensity","aggregationPropensityInVivo","hydrophobicity","hydrophobicMoment","isoelectricPoint"]
d=[]
newds={"id":list(dic.keys())}
fMinMax=pd.read_csv("ampNonAMPminMaxesTable50.csv")
for j in features:
    p=[]
    for i in p6.loc[:,j]:
        p.append((i-fMinMax.loc[0, j])/(fMinMax.loc[1, j]-fMinMax.loc[0, j]))
    newds[j]=p 
newdf=pd.DataFrame(data=newds)
newdf
rfModel = pickle.load(open('rf_AMP_classifier_final_model.sav', 'rb'))
#svcModel = pickle.load(open('svc_AMP_classifier_final_model.sav', 'rb'))
mr=newdf.iloc[:, 0].values
print('-------AMP/non-AMP Classification Results-------')
print('Random Forest Classifier')
staphActCheckList=[]
staphActCheckListNo=[]
for i in mr:
    take=np.where(mr == i)
    pred=rfModel.predict(newdf.iloc[take[0],1:].values)
    #print(i,pred[0])
    print(i,labelDict[pred[0]],dsss['seq'][i-1])
    if labelDict[pred[0]]=='AMP':
        staphActCheckList.append(dsss['seq'][i-1])
        staphActCheckListNo.append(i)




#from __future__ import print_function
import datetime
from isoelectric import ipc #Requires Isoelectric Package
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import csv
import statistics
from propy import PyPro #Requires Installing Propy Package
import propy #Requires Installing Propy Package
from propy.PyPro import GetProDes #Requires Installing Propy Package
import pickle
labelDict={"Neg":"Negative","Plus":"Positive"}
#Defining Normalized Hydrophobicity Calculators functions

import argparse
import math
import os
import time


aa_charge = {'E': -1, 'D': -1, 'K': 1, 'R': 1}

#

def Merge(dict1, dict2): 
    return(dict2.update(dict1))
pp=[]

#Create Dictionary of Sequences
c={}
x = range(1,len(staphActCheckList)+1)
dic={}
for i in staphActCheckList:
    dic[staphActCheckList.index(i)+1]=i
#Calculate Properties
aggDict = pickle.load(open('aggDict.pkl', 'rb'))#Requires aggDict File
chargeDict = pickle.load(open('chargeDict.pkl', 'rb'))#Requires chargeDict File
#isoElP=[]
netcharge=[]
agg=[]
chargeDensities=[]
mw=[]
#normalizedHydrophobicity=[]
#normalizedHydrophobicMoment=[]
dsssb={"seq":staphActCheckList}
for key, val in dic.items():
    #isoElP.append(ipc.predict_isoelectric_point(val,"IPC_peptide"))
    agList=[]
    chargeList=[]
    mwList=[]
    for j in list(val):
        agList.append(aggDict[j])
        chargeList.append(chargeDict[j])
    agg.append(statistics.mean(agList))
    chargeDensities.append(sum(chargeList)/ipc.calculate_molecular_weight(i))
    mw.append(ipc.calculate_molecular_weight(val))
    netcharge.append(sum(chargeList))
    #normalizedHydrophobicity.append(round(sum(assign_hydrophobicity(val))*(-1)/len(val),2))
    #normalizedHydrophobicMoment.append(round(calculate_moment(assign_hydrophobicity(val)),2))
    Des=GetProDes(val)
    aACompDes=Des.GetAAComp() #AA Composition
    dPCompDes=Des.GetDPComp() #DP Composition
    moreauBrotoAutoDes=Des.GetMoreauBrotoAuto() #
    moranAutoDes=Des.GetMoranAuto() #
    gearyAutocDes=Des.GetGearyAuto() #
    cTDDes=Des.GetCTD()
    aPAACDes=Des.GetAPAAC(lamda=5,weight=0.05) #Type II PseudoAAC
    paacDes=Des.GetPAAC(lamda=5,weight=0.05) #Type I PseudoAAC
    sOCNDes=Des.GetSOCN()
    qSODes=Des.GetQSO(maxlag=30,weight=0.1)
    Merge(qSODes, sOCNDes)
    Merge(sOCNDes, paacDes)
    Merge(paacDes, aPAACDes)
    Merge(aPAACDes, cTDDes)
    Merge(cTDDes, gearyAutocDes)
    Merge(gearyAutocDes, moranAutoDes)
    Merge(moranAutoDes, moreauBrotoAutoDes)
    Merge(moreauBrotoAutoDes, dPCompDes)
    Merge(dPCompDes, aACompDes)
    pddd=aACompDes.copy()
    pp.append(pddd)
dsssb["molecularWeight"]=mw
dsssb["netCharge"]=netcharge
dsssb["chargeDensity"]=chargeDensities
dsssb["AggregationPropensityInVivo"]=agg
#dsss["aggregationPropensityInVivo"]=agg
#dsss["chargeDensity"]=chargeDensities
df=pd.DataFrame(data=pp)
p5=pd.DataFrame(data=dsssb)
p6=pd.concat([p5, df.reindex(p5.index)], axis=1)
#Rescaling Data
#Needs Final Features MinMaxes
p=[]
features=['netCharge','molecularWeight','chargeDensity','AggregationPropensityInVivo','MoreauBrotoAuto_ResidueVol1','MoreauBrotoAuto_ResidueVol5','MoreauBrotoAuto_Mutability1','MoreauBrotoAuto_Mutability4','MoreauBrotoAuto_Mutability8','MoreauBrotoAuto_Mutability12','MoranAuto_FreeEnergy7','MoranAuto_FreeEnergy11','MoranAuto_FreeEnergy21','MoranAuto_ResidueVol11','MoranAuto_Steric4','MoranAuto_Steric8','MoranAuto_Steric9','MoranAuto_Steric11','MoranAuto_Mutability11','MoranAuto_Mutability22','GearyAuto_Polarizability4','GearyAuto_FreeEnergy19','GearyAuto_Steric4','GearyAuto_Steric6','GearyAuto_Steric9','_SolventAccessibilityC1','_SolventAccessibilityC3','_ChargeC3','_PolarityC1','_PolarizabilityT13','_PolarizabilityD3001','_SecondaryStrD2001','_ChargeD1075','_PolarityD1001','PAAC22','PAAC25','PAAC26','PAAC29','tausw1','tausw2','tausw3','tausw4','tausw5','tausw8','tausw9','tausw13','QSOSW17','QSOSW24','QSOSW25','QSOSW29','QSOSW33']
d=[]
newds={"id":list(dic.keys())}
fMinMax=pd.read_csv("staphDsFeaturesMinMax.csv")
for j in features:
    p=[]
    for i in p6.loc[:,j]:
        p.append((i-fMinMax.loc[0, j])/(fMinMax.loc[1, j]-fMinMax.loc[0, j]))
    newds[j]=p 
newdf=pd.DataFrame(data=newds)
newdf
#rfModel = pickle.load(open('rf_AfterFS2-51f_final_model.sav', 'rb'))
eclf5Model = pickle.load(open('eclf5.sav', 'rb'))
mr=newdf.iloc[:, 0].values
print('-------S.Aureus Activity Prediction-------')
print('Hybrid Classifier Results')
for i in mr:
    take=np.where(mr == i)
    pred=eclf5Model.predict(newdf.iloc[take[0],1:].values)
    print(staphActCheckListNo[i-1],labelDict[pred[0]],staphActCheckList[i-1])