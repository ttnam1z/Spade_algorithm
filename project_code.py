#!curl -L "http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/diab_trans.data" > diab_trans.csv
import numpy as np
import pandas as pd
import scipy as sp
import os

#CDir = os.getcwd()
#diab_data=pd.read_csv(CDir + "\\Metro_Interstate_Traffic_Volume.csv")
#print(diab_data.columns)
#print(diab_data.values[diab_data['code']=="id_65",:])
#print(diab_data.code[4:7])

def Spade(inputItem, input_data){
    #convert to number all data

    minSup=1
    #1.Find frequent items
    #Count support for items
    for each sid
        for each eid
            get list items in eid
            if not in matched list
                item count++
                put item in matched list 
        clear matched list
    for each item in item list
        if count >minSup
            add to frequent items list
        else 
            remove items in eid in sid
    #2.Find frequent 2-sequences
    for each sid
        for each item in sid
            create 2-sequence
            if not in matched 2-seq list
                2-sequence count++
                put 2-sequence in 2-seq list
        clear matched 2-seq list

    for each item in 2-seq list
        if count >minSup
            add to frequent 2-seq items list    


    #Enumerating frequent sequences
    Enumerating_frequent(frequentlist)
}

def Enumerating_frequent(frequentlist){
    for Ai in frequentlist
        set Ti list empty
        for Aj in frequentlist j>i
            R = Ai join Aj
            if(Prune(R) == FALSE)
                id-list(R) = id-list(Ai) join id-list(Aj)
                if count of R > minsup
                    Add R to Ti list
                    

        if Depth-first-search
            Enumerating_frequent(Ti)
    if Breadth-first-search
        for all Ti != empty
            Enumerating_frequent(Ti)      
}

def Prune(seq){
    for sub sequence s of seq
        if s has processed and s not belong to frequent list
            return True
    return False
}

def preprocess(data){

}
