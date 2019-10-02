# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:50:20 2019

@author: Barry Gan
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 22:56:20 2019

@author: Barry Gan
"""

import numpy as np
import csv as csv
import matplotlib.pyplot as plt
import pandas as pd
import statistics as stc

#read in csv file as array
##the read in funtion is cited from https://gist.github.com/dettmering/3767366
def readcsv(filename):
    ifile = open(filename, "r")
    reader = csv.reader(ifile, delimiter=",")

    rownum = 0
    a = []

    for row in reader:
        a.append(row)
        rownum += 1

    ifile.close()
    return a

#define a class, in which sample and corresponding data is integrated
class sample_data:
    def __init__(self,samplelist,conditionlist:list):
        self.barcode=''
        self.data=[]
        self.data_name=[]
        self.name=''
        self.theo_rep=0
        self.prac_rep=0
        self.sample_list=samplelist
        self.condition_list=conditionlist
        self.namecode=''
        return
    def write_name(self,barcode:str):
        self.barcode=barcode
        return
    def write_data(self,name,data):
        self.data.append(data)
        self.data_name.append(name)
        return
    def decode(self):
        #get the name, which theo and prac rep of the sample is
        #spit the barcode

        sample_list=self.sample_list
        if self.barcode=='':
            self.name='null'
        else:
            name_string=self.barcode.split("-")
            self.theo_rep=int(name_string[1])
            self.prac_rep = int(name_string[2])
            a=name_string[0]
            self.namecode=a
            condition_no=len(self.condition_list)
            for i in range(0,condition_no):
                if a==sample_list[i]:
                    self.name=self.condition_list[i]
        return
#define a class, in which all samples (rep) of 1 certain condition are clustered
#class condition_cluster:
## theo_rep, prac_rep is the rep considered without outlier
class condition_cluster:
    def __init__(self,exp_name,namecode,name,theo_rep:int,prac_rep:int):
        self.theo_rep=theo_rep
        self.prac_rep=prac_rep
        self.exp_name=exp_name
        self.namecode=namecode
        self.name=name
        self.samples=[ [ None for y in range(self.prac_rep) ] for x in range(self.theo_rep)]

        self.data570=[ [ None for y in range(self.prac_rep)] for x in range(self.theo_rep)]
        self.data600=[ [ None for y in range(self.prac_rep)] for x in range(self.theo_rep)]

#the parameters, here specified for Ablue
        self.oxi570=80586
        self.oxi600=117216
        self.red570=155677
        self.red600=14652

#initialize variable to store value, here specified for Ablue
        self.redablue=[]
        self.cdtavr570=[]
        self.cdtavr600=[]
        self.redavr=0
        self.redsd=0
        return
    #specific for Ablue cal the average of prac rep
    def prac_rep_avr(self):
        print(self.name)
        for i in range(0,self.theo_rep):
            self.cdtavr570.append(np.mean(self.data570[i]))
            self.cdtavr600.append(np.mean(self.data600[i]))
            #print(self.data570[i])
            #self.cdtavr600.append(np.mean(self.data600[i]))
        return
    #cal reduction of alamar blue reagent
    #specified for Ablue
    def red_ablue_rate(self,c600,c570):
        #test
        print(self.name)
        eoxi600=self.oxi600
        eoxi570=self.oxi570
        ered600=self.red600
        ered570=self.red570
        theo_rep=self.theo_rep
        for i in range(0,theo_rep):
            #test
            print("theo_rep "+str(i+1))
            a570=self.cdtavr570[i]
            a600=self.cdtavr600[i]
            redablue=((eoxi600*a570)-(eoxi570*a600))/((ered570*c600)-(ered600*c570))
            print('%.2f%%' % (redablue * 100))
            self.redablue.append(redablue)
        self.redavr=np.mean(self.redablue)
        self.redsd=np.std(self.redablue)
        return
        

    



#this function log the name of each sample,from the well-plan csv
def write_name(samplelist,conditionlist:list,well_plan_csv:str):
    well_plan = readcsv(well_plan_csv)
    sample_set=[]
    for i in range(0,8):
        for j in range(0,12):
            one_data=sample_data(sample_list,conditionlist)
            one_data.write_name(well_plan[i][j])
            one_data.decode()
            sample_set.append(one_data)
    return sample_set

#this function write datas under each sample, from the output of the csv samples
def write_data(sample_set,var_name,data_csv:str):
    data_frame = readcsv(data_csv)
    for i in range(0,8):
        for j in range(0,12):
            no=i*12+j
            sample_set[no].write_data(var_name,float(data_frame[i][j]))
    return sample_set
#the corresponding condition code and the codition
sample_list=['0','1','2','3','4','5','6','7']
conditionlist = ['- ctrl','c1+','c1-','c2+','c2-','Col','scfd +','scfd -']
experiment_name='5d'


sample_set=write_name(sample_list,conditionlist,'well_plan5d.csv')
sample_set=write_data(sample_set,'5d570nm','data5d570.csv')
sample_set=write_data(sample_set,'5d600nm','data5d600.csv')

#cluster samples according to each condition and collect in a list
def Alamar_blue_cluster(exp_name,sample_list,condition_list,sample_set,theo_rep,prac_rep):
    a=[]
    length=len(sample_list)
#initialize conditions, into which sample data would be fed in
#initialize the - ctrl
    for i in range(0,1):
        b=condition_cluster(exp_name,sample_list[i],condition_list[i],2,2)
        a.append(b)

#initialize all other conditions
    for i in range(1,length):
        b=condition_cluster(exp_name,sample_list[i],condition_list[i],3,2)
        a.append(b)
 
#cluster the data according to the conditions
    for i in sample_set:
        #get rid of null well
        if i.name!='null':
            namecode=int(i.namecode)
            theo_rep=i.theo_rep
            prac_rep=i.prac_rep
            if theo_rep<=a[namecode].theo_rep and prac_rep<=a[namecode].prac_rep:
                (a[namecode].samples)[theo_rep-1][prac_rep-1]=i
                (a[namecode].data570)[theo_rep-1][prac_rep-1]=i.data[0]
                (a[namecode].data600)[theo_rep-1][prac_rep-1]=i.data[1]

                
            else:
                continue
    return a
Ablue5d_cdt=Alamar_blue_cluster(experiment_name,sample_list,conditionlist,sample_set,3,2);

#for i in range(0,8):

#from here, specified for Ablue
#get c570 c600
c600=np.mean([Ablue5d_cdt[0].data600[0][0],Ablue5d_cdt[0].data600[0][1],Ablue5d_cdt[0].data600[1][0],Ablue5d_cdt[0].data600[1][1]])
c570=np.mean([Ablue5d_cdt[0].data570[0][0],Ablue5d_cdt[0].data570[0][1],Ablue5d_cdt[0].data570[1][0],Ablue5d_cdt[0].data570[1][1]])
length=len(sample_list)
print("c600",c600)
print("c500",c570)

#calculate reduction rate of alamar blue
#write the result in x,y,sd(for error bar) for later plotting of the result
y=[]
x=conditionlist[1:8]
sd=[]
for i in range(0,length):
    Ablue5d_cdt[i].prac_rep_avr()
    Ablue5d_cdt[i].red_ablue_rate(c600,c570)
    if i>0:
        y.append(Ablue5d_cdt[i].redavr*100)
        sd.append(Ablue5d_cdt[i].redsd*100)

#plot, cited from https://blog.csdn.net/songyunli1111/article/details/83625639
index = np.arange(7)
values = y
SD = sd
plt.title('5d Reduction of Alamar Blue Reagent')
plt.ylabel('Reduction rate of Alamar Blue/%')
plt.bar(index, values, yerr = SD, error_kw = {'ecolor' : '0.2', 'capsize' :6}, alpha=0.7)
plt.xticks(index+0.2,x)
plt.legend(loc=2)
plt.savefig("result5d.svg")
plt.show()
plt.clf()
plt.cla()
plt.close()
#output result in csv for manually check
header=['Condition','Theo_rep','570nm absorb of Prac_rep','Avr 570 absor','600nm absorb of Prac_rep','Avr 600 absor','Ablue Redc Rate','Avr Ablue Redc Rate','c570','c600']
lines=list()
for i in range(1,length):
    theo_rep=Ablue5d_cdt[i].theo_rep
    for j in range(0,3):
        line=list()

        line.append([Ablue5d_cdt[i].name])
        line.append(str(j+1))
        line.append(Ablue5d_cdt[i].data570[j])
        line.append(Ablue5d_cdt[i].cdtavr570[j])
        line.append(Ablue5d_cdt[i].data600[j])
        line.append(Ablue5d_cdt[i].cdtavr600[j])
        line.append(Ablue5d_cdt[i].redablue[j])
        line.append(Ablue5d_cdt[i].redavr)
        line.append(c570)
        line.append(c600)
        lines.append(line)
with open("result5d.csv", "w", newline='') as f:
    writer = csv.writer(f, delimiter=',')
    writer.writerow(header) # write the header
    # write the actual content line by line
    for l in lines:
        writer.writerow(l)
