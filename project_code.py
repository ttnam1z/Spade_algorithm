#!curl -L "http://staff.ii.pw.edu.pl/~gprotazi/dydaktyka/dane/diab_trans.data" > diab_trans.csv
import numpy as np
import pandas as pd
import scipy as sp
import os

#sequenclist structure:
#   items: list of item.
#   num: number of item in list
# 
#item structure: 
#   name: list of event. each event is list of item
#   pairInfo: matrix has 2 column : sid,eid
#
class SequenceList:
    def __init__(self):
        self.items=[]
        self.num=0
    def add(self,item,sub=1):
        #search list
        fi = False
        for idx in range(self.num):
            length=len(item["name"])
            if(len(self.items[idx]["name"]) == length):
                match= True
                for cnt in range(length):
                    if (set([tuple(itemset) for itemset in self.items[idx]["name"][cnt]])) \
                        != (set([tuple(itemset1) for itemset1 in item["name"][cnt]])):
                        match = False
                        break
                if match:  # found sequence in list, add pair Info
                    fi = True
                    self.items[idx]["pairInfo"]=np.vstack([self.items[idx]["pairInfo"],item["pairInfo"]])
                    break
        if not fi: # add new sequence
            pairInfo = np.empty([0,2])
            pairInfo = np.vstack([pairInfo,item["pairInfo"]])
            new_item = {"name":item["name"],"pairInfo":pairInfo, "sub":sub}
            self.items.append(new_item)
            self.num += 1
    def remove(self,item):
        self.items.remove(item)
        self.num = len(self.items)
    def getInfo(self,event):
        pairInfo = np.empty([0,2])
        for idx in range(self.num):
            length=len(event)
            if(len(self.items[idx]["name"]) == length):
                match= True
                for cnt in range(length):
                    if (set([tuple(itemset) for itemset in self.items[idx]["name"][cnt]])) \
                        != (set([tuple(itemset) for itemset in event[cnt]])):
                        match = False
                        break
                if match:  # found sequence in list, add pair Info
                    pairInfo=np.vstack([pairInfo,self.items[idx]["pairInfo"]])
                    break
        return pairInfo

#SidList structure:
#   seqs: list of sequence
#   num: number of seqs
#Sequence structure:
#   id: sequence id
#   pairInfo:  
#       name: list of events. event is list of items
#       eid: event id
#
class SidList:
    def __init__(self):
        self.seqs=[]
        self.num=0
    def add(self,seq):
        #search list
        fi = False
        for idx in range(self.num):
            if self.seqs[idx]["id"]== seq["id"]:
                fi = True
                self.seqs[idx]["pairInfo"].append({"item":seq["name"], "eid":seq["eid"]})
                break
        if not fi:
            new_seq = {"id":seq["id"],"pairInfo":[{"item":seq["name"], "eid":seq["eid"]}]}
            self.seqs.append(new_seq)
            self.num += 1

    def remove(self,seq):
        self.seqs.remove(seq)
        self.num = len(self.seqs)

def SeparateStr(str):
    #a,b->c->d  : ['a,b->c', 'd']
    #a->b->c,d   : ['a->b->c', 'd']
    comma = str.rfind(',')
    imply = str.rfind('->')
    if comma > imply:
        return str.rsplit(',',1)
    else:
        return str.rsplit('->',1)

def Compare2Seqs(seq1,seq2): 
    #seq is list of events. each event is list of items.
    #True is match, False is unmatch

    match = True
    length=len(seq2)
    if(len(seq1) == length):
        for cnt in range(length):
            if (set([tuple(itemset) for itemset in seq1[cnt]])) \
                != (set([tuple(itemset) for itemset in seq2[cnt]])):
                match = False
                break
    else:
        match = False
    return match

def FindPairInfo(infoLs1, infoLs2, tp1=-1):
    #Find mix pair info of 2 id list info
    mixInfo = np.empty([0,2])
    sid1 = infoLs1.shape[0]
    sid2 = infoLs2.shape[0]
    if tp1 < 0:
        for info1 in range(sid1):
            for info2 in range(sid2):
                if infoLs1[info1,0] == infoLs2[info2,0] and infoLs1[info1,1] < infoLs2[info2,1]:
                    mixInfo = np.vstack([mixInfo,infoLs2[info2,:]])
    elif  tp1 > 0:
        for info1 in range(sid1):
            for info2 in range(sid2):
                if infoLs1[info1,0] == infoLs2[info2,0] and infoLs1[info1,1] > infoLs2[info2,1]:
                    mixInfo = np.vstack([mixInfo,infoLs1[info1,:]])
    else:
        for info1 in range(sid1):
            for info2 in range(sid2):
                if infoLs1[info1,0] == infoLs2[info2,0] and infoLs1[info1,1] == infoLs2[info2,1]:
                    mixInfo = np.vstack([mixInfo,infoLs1[info1,:]])
    return mixInfo

def copyLs(inLs):
    #use for list of list to avoid reference 
    ls = []
    for item in inLs:
        new_item = item.copy()
        ls.append(new_item)
    return ls

def Join2Seqs(seq1, seq2):
    #list joined sequences
    #seq1 = {"name":[["b","a"]]}
    #seq2 = {"name":[["a","c"]]}
    
    itemLs = []

    #join 2 items
    l1 = len(seq1["name"])
    l2 = len(seq2["name"])

    if l2 < (l1 - 1) or l2 > (l1 + 1):
        return itemLs

    if l1 == 1:
        if l2 == 1:
            inter_ = set(seq1["name"][-1]).intersection(set(seq2["name"][-1]))
            dif_1 = set(seq1["name"][-1]).difference(set(seq2["name"][-1]))
            dif_2 = list(set(seq2["name"][-1]).difference(set(seq1["name"][-1])))
            if len(dif_1) == len(dif_2):
                if len(dif_1) == 1:
                    # ac join ab
                    #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],0) 
                    name = copyLs(seq1["name"])
                    name[-1].append(dif_2[0])
                    #item = {"name":name,"pairInfo":pairInfo}
                    itemLs.append({"name":name})
                elif  len(dif_1) == 0:
                    #abc join abc
                    for idx in range(len(seq1["name"][-1])):
                        #pairInfo = FindPairInfo(seq1["pairInfo"],seq1["pairInfo"],-1) 
                        name = copyLs(seq1["name"])
                        name.append([seq1["name"][-1][idx]])
                        #item = {"name":name,"pairInfo":pairInfo}
                        itemLs.append({"name":name})
        elif  l2 == 2:
            inter_ = set(seq1["name"][-1]).intersection(set(seq2["name"][-2]))
            dif_1 = set(seq1["name"][-1]).difference(set(seq2["name"][-2]))
            if len(inter_) == len (seq2["name"][-2]) and len(dif_1) == 1 and len(seq2["name"][-1]) == 1:
                #joinable
                #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],-1) 
                name = copyLs(seq1["name"])
                name.append([seq2["name"][-1][0]])
                #item = {"name":name,"pairInfo":pairInfo}
                itemLs.append({"name":name})
    else:
        if l1 == l2:
            inter_ = set(seq1["name"][-1]).intersection(set(seq2["name"][-1]))
            dif_1 = list(set(seq1["name"][-1]).difference(set(seq2["name"][-1])))
            dif_2 = list(set(seq2["name"][-1]).difference(set(seq1["name"][-1])))
            if len(dif_1) == len(dif_2) and Compare2Seqs(seq1["name"][0:-1],seq2["name"][0:-1]):
                if len(dif_1) == 1: 
                    #joinable
                    #
                    #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],0) 
                    name = copyLs(seq1["name"])
                    name[-1].append(dif_2[0])
                    #item = {"name":name,"pairInfo":pairInfo}
                    itemLs.append({"name":name})
                    #
                    #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],-1) 
                    name = copyLs(seq1["name"])
                    name.append([dif_2[0]])
                    #item = {"name":name,"pairInfo":pairInfo}
                    itemLs.append({"name":name})
                    #
                    #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],1) 
                    name = copyLs(seq2["name"])
                    name.append([dif_1[0]])
                    #item = {"name":name,"pairInfo":pairInfo}
                    itemLs.append({"name":name})
                elif  len(dif_1) == 0:
                    #abc join abc
                    for idx in range(len(seq1["name"][-1])):
                        #pairInfo = FindPairInfo(seq1["pairInfo"],seq1["pairInfo"],-1) 
                        name = copyLs(seq1["name"])
                        name.append([seq1["name"][-1][idx]])
                        #item = {"name":name,"pairInfo":pairInfo}
                        itemLs.append({"name":name})
                
        elif  l1 < l2:
            inter_ = set(seq1["name"][-1]).intersection(set(seq2["name"][-2]))
            dif_1 = set(seq1["name"][-1]).difference(set(seq2["name"][-2]))
            if len(inter_) == len (seq2["name"][-2]) and len(dif_1) == 1 and len(seq2["name"][-1]) == 1 \
                and Compare2Seqs(seq1["name"][0:-1],seq2["name"][0:-2]):
                #joinable
                #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],-1)  
                name = copyLs(seq1["name"])
                name.append([seq2["name"][-1][0]])
                #item = {"name":name,"pairInfo":pairInfo}
                itemLs.append({"name":name})
        else: #l1>l2
            inter_ = set(seq2["name"][-1]).intersection(set(seq1["name"][-2]))
            dif_1 = set(seq2["name"][-1]).difference(set(seq1["name"][-2]))
            if len(inter_) == len (seq1["name"][-2]) and len(dif_1) == 1 and len(seq1["name"][-1]) == 1\
                and Compare2Seqs(seq2["name"][0:-1],seq1["name"][0:-2]):
                #joinable
                #pairInfo = FindPairInfo(seq1["pairInfo"],seq2["pairInfo"],1) 
                name = copyLs(seq2["name"])
                name.append([seq1["name"][-1][0]])
                #item = {"name":name,"pairInfo":pairInfo}
                itemLs.append({"name":name})
    return itemLs

def JoinFreqItem(eid1, eid2, sid):
    # eid of eid2 need greater or equal eid of eid1
    if eid2["eid"] < eid1["eid"]:
        return []
    ls = []
    #mix 2 items
    if eid1["item"] == eid2["item"]:
        ls = [{"name":[eid1["item"][0],eid2["item"][0]],"pairInfo":[sid,eid2["eid"]]}]
    else:
        if eid2["eid"] == eid1["eid"]:
            ls = [{"name":[[eid1["item"][0][0], eid2["item"][0][0]]],"pairInfo":[sid,eid2["eid"]]}]
        else:
            ls = [{"name":[eid1["item"][0],eid2["item"][0]],"pairInfo":[sid,eid2["eid"]]}]
    return ls

def FindEquiCls(seqSet, freqItems):
    #only find equivalance class for k = 1
    k = 1
    equiClsLs = []
    for idx in range(freqItems.num):
        equiClsLs.append([])
        for jdx in range(seqSet.num):
            inters = set(freqItems.items[idx]["name"][0]).intersection(set(seqSet.items[jdx]["name"][0]))
            dif1 = set(freqItems.items[idx]["name"][0]).difference(inters)
            if len(inters) == k and len(dif1) == 0: 
                #add list
                equiClsLs[idx].append(seqSet.items[jdx])
    return equiClsLs

def Enumerating_frequent(equiCls, F, prs, method, minSup, in_idx, maxLen=5):
    T = []
    cur_idx = in_idx + 1
    if cur_idx >= 5:
        return
    if len(F) == cur_idx:
        F.append(SequenceList())
        prs.append([])
    num = len(equiCls) 
    for idx in range(num):
        T.append([])
    for idx in range(num):
        for jdx in range(idx,num):    
            #join 2 seq
            ls = Join2Seqs(equiCls[idx],equiCls[jdx])
            #check frequent sequence
            for item in ls:
                if (not CheckProcessed(item,prs[cur_idx-2])):
                    if Prune(item, F, cur_idx-1)==False:
                        #count sup
                        pInfo1 = np.empty([0,2])
                        for event in item["name"]:
                            len1 = len(event)
                            if len1 == cur_idx+1: # fix cant find when event is not listed in frequence list. example: seq = [[abc]]
                                pInfo2 = FindPairInfo(F[len1-2].getInfo([event[0:-1]]),F[0].getInfo([[event[-1]]]),0) #fix missing some special seq
                            else:
                                pInfo2 = F[len1-1].getInfo([event])
                            if len(pInfo1) == 0:
                                pInfo1 = pInfo2
                            else:
                                pInfo1 = FindPairInfo(pInfo1, pInfo2)
                        sup= np.unique(pInfo1[:,0]).shape[0]
                        if sup > minSup:
                            T[idx].append(item)
                            if idx != jdx:
                                T[jdx].append(item)
                            seq = {"name":item["name"],"pairInfo":pInfo1}
                            F[cur_idx].add(seq,sup)
                    prs[cur_idx-2].append(item)
        if method == "depth":
            if len(T[idx]) > 0:
                Enumerating_frequent(T[idx],F,prs,method,minSup,cur_idx)
    if method == "breadth":
        for item in T:
            if len(item) > 0:
                Enumerating_frequent(item,F,prs,method,minSup,cur_idx)
    

def GetSubsequence(seq):
    #seq = {"name":[["b"]]}
    #seq = {"name":[["c","d","e","f"],["b"]]}
    #seq = {"name":[["c"],["b"]]}
    #seq = {"name":[["a","b"],["c","d","e","f"],["b"]]}
    #seq = {"name":[["a","b"],["c"],["b"]]}
    #seq = {"name":[["a","b"],["c","d","e","f"],["b","a"]]}
    #seq = {"name":[["c","d","e","f"]]}
    ls=[]
    len1 = len(seq["name"])
    len2 = len(seq["name"][-1])
    prefix = seq["name"][0:-2]
    if len2 == 1:
        if len1 == 1:
            ls.append(seq["name"])
        elif  len1 == 2:
            ls.append(seq["name"][0:-1])
            len3 = len(seq["name"][-2])
            for idx in range(len3):
                item = seq["name"][-2].copy()
                item.pop(idx)
                if len(item) > 0:
                    ls.append([item,seq["name"][-1]])
                else:
                    ls.append([seq["name"][-1]])
        elif  len1 > 2:
            ls.append(seq["name"][0:-1])
            len3 = len(seq["name"][-2])
            for idx in range(len3):
                item = seq["name"][-2].copy()
                item.pop(idx)
                seq1 = seq["name"][0:-2]
                if len(item) > 0:
                    seq1.append(item)
                seq1.append(seq["name"][-1])
                ls.append(seq1)
    elif  len2 > 1: 
        if len1 == 1:
            for idx in range(len2):
                item = seq["name"][-1].copy()
                item.pop(idx)
                ls.append([item])
        elif  len1 > 1:
            for idx in range(len2):
                item = seq["name"][-1].copy()
                item.pop(idx)
                seq1 = seq["name"][0:-1]
                seq1.append(item)
                ls.append(seq1)
    return ls

def CheckProcessed(seq,fSet):
    find = False
    for idx in range(len(fSet)):
        if Compare2Seqs(seq["name"],fSet[idx]["name"]):
            find = True
            break
    return find

def Prune(seq, F, index):
    #check subsequence of seq is a frequent sequence
    find = 0
    ls = GetSubsequence(seq)
    for sSeq in ls:
        for idx in range(F[index].num):
            if Compare2Seqs(sSeq,F[index].items[idx]["name"]):
                find +=1
                break
    return (find != len(ls))

#data structure
#item
#(sid, eid)
def Spade(inSup, input_data, method):
    #convert to number all data
    F=[]
    l1 = np.unique(input_data[:,0])
    minSup= int(inSup * np.unique(input_data[:,0]).shape[0])
    #minSup = 0
    #1.Find frequent items
    #create vertical data
    itemLs = SequenceList()
    le = input_data.shape[0]
    for idx in range(le):
        sid = int(input_data[idx,0])
        eid = int(input_data[idx,1])
        items = input_data[idx,2].split(" ")
        
        for it in items:
            item = {"name":[[it]],"pairInfo":[sid,eid]}
            itemLs.add(item)

    #Count support for items
    #s_num=np.unique(input_data[0]).shape[0]
    F.append(SequenceList())
    for item in itemLs.items:
        item_sup = np.unique(item["pairInfo"][:,0]).shape[0]
        if item_sup > minSup:
            F[0].add(item,item_sup)

    #2.Find frequent 2-sequences
    #Convert to horizon data
    sidLs = SidList()
    for item in F[0].items:
        for idx in range(item["pairInfo"].shape[0]):
            seq = {"name":item["name"], "id":item["pairInfo"][idx,0], "eid":item["pairInfo"][idx,1]}
            sidLs.add(seq)
    
    #Create list 2-sequences
    item2Ls = SequenceList()
    for seq in sidLs.seqs:
        pair_len = len(seq["pairInfo"])
        for idx in range(pair_len):
            for jdx in range(idx+1,pair_len):
                #mix 2 items
                if (seq["pairInfo"][jdx]["eid"] >= seq["pairInfo"][idx]["eid"]):
                    for item in JoinFreqItem(seq["pairInfo"][idx],seq["pairInfo"][jdx],seq["id"]):
                        item2Ls.add(item)
                else:
                    for item in JoinFreqItem(seq["pairInfo"][jdx],seq["pairInfo"][idx],seq["id"]):
                        item2Ls.add(item)
    #Count support
    F.append(SequenceList())
    for item in item2Ls.items:
        item_sup = np.unique(item["pairInfo"][:,0]).shape[0]
        if item_sup > minSup:
            F[1].add(item,item_sup)
    
    prs = [] # check processed
    #Enumerating frequent sequences
    for equiCls in FindEquiCls(F[1],F[0]):
        if len(equiCls) > 0:
            Enumerating_frequent(equiCls, F, prs, method, minSup, 1)
    return F


def preprocess(data):
    #TODO
    pass

def main():
    CDir = os.getcwd()
    data = np.empty([0,3])
    k=0
    #with open(CDir + "\\tags_full.data",encoding="utf8",mode="r") as f:
    with open(CDir + "\\test.data",encoding="utf8",mode="r") as f:
        for line in f:
            ls = line.strip("\n").split(" ", 3)
            data = np.vstack([data,[ls[0],ls[1],ls[3]]])
            k +=1
            #if k>100:
            #    break

    preprocess(data)
    F = Spade(0.3,data,"depth")
    for item in F:
        print(item.num)
        print("\n")
main()
