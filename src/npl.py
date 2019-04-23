#!/usr/bin/python 
#Author: Linhai Zhao
#NPL scoring functions to detect over-sharing in affecteds in extended families
from __future__ import division
#from memory_profiler import profile
import argparse
import os
import sys
import re
import random
import copy
from zipfile import ZipFile
from ctypes import *
from multiprocessing import Process, Value, Manager, Pool
from multiprocessing.sharedctypes import Array
import multiprocessing
from itertools import repeat,product,permutations
import math
import numpy as np
import time
import string
from scipy import stats
from scipy.optimize import fsolve
from sympy.solvers import solve,solveset
from sympy import Symbol,Interval
from RVNPLcpp import cmissingparents,cmissing_infer,sall_cpp,cpostInv
from RVNPL import screen_output,ibd_module
#################
manager=Manager()

def random_dist(dic):
    r = random.uniform(0,1)
    s = 0
    for item,prob in dic.iteritems():
        s += prob
        if s >= r:
            return item
    return item

def ransample(args):
        structnum,fam_struct = args
        t_ibd = [0,0]
        fs_num = len(fam_struct)
        for i in range(fs_num):
            s_ibd = [random_dist(fam_struct[i]) for tmp in xrange(structnum[i])]
            t_ibd[0] += sum([x[0] for x in s_ibd])
            try:
                t_ibd[1] += sum([x[1] for x in s_ibd])
            except TypeError:
                t_ibd[1] = None
        t_ibd[0]=float('%.9f'%(t_ibd[0]))
        try:
            t_ibd[1]=float('%.9f'%(t_ibd[1]))
        except:
            pass
        return tuple(t_ibd)

def null_inv(args):
#calculate null ibd by given inheritance vector
    v,fam=args
    tmp = bin(v)[2:]                     #get binary format
    inv = '0'*(2*len(fam.invnf)-len(tmp))+tmp
    inv = [int(tmpv) for tmpv in inv]
    fam_alleles=fam.assign_allele(inv)
    ibd_n=0 #sum of pairwise IBD among affected pairs
    for aff_pair in fam.affected_pairs:
        ibd_pair=len(set(fam_alleles[aff_pair[0]])&set(fam_alleles[aff_pair[1]]))
        ibd_n+=ibd_pair
    tmp_sall=fam.sall(inv,fam_alleles,[])
    return ibd_n,tmp_sall,v


def null_generator(procID,founderid,fam,queue,all_ibd,null_ibd,pinv_sall,pinv_key_dict,pinv_pair_dict,sall_flag,infer_flag):
    outID = 0
    while True:
        try:
            data = queue.get()
            if data is None:
                break
            else:
                ibd_result=[None,None]
                foundergt_list = []
                foundergt = {}
                new_alleles = {}
                nfnum = len(fam.nonfounder)
                if isinstance(data,tuple):
                    #calculate null ibd by genotypes
                    v,fgt=data
                    for idx,f in enumerate(fam.founder):
                        foundergt[f] = list(fgt[2*idx:2*idx+2])
		else:
                    #permutations
                    random_gt=None
                    if fam.missing_all != []:
                        tmp_dic=fam.conditional_prob['~combined'][0]
                        gt_id=fam.conditional_prob['~combined'][1]
                        r = random.uniform(0,1)
                        s = 0
                        for gt,prob in tmp_dic.iteritems():
                            s += prob
                            if s >= r:
                                random_gt=gt
                                break
                        if random_gt == None:
                            random_gt = tmp_dic.keys()[0]
                    for f in founderid:
                        if f in fam.missing_all:
                            idx=gt_id.index(f)
                            foundergt[f]=[ord(x)-96 for x in random_gt[2*idx:2*idx+2]]
                        else:
                            foundergt[f]=fam.fam_dict[f]['gt'][2*fam.mid:2*fam.mid+2]
                    #Randomly assign an inheritance vector
                    v = random.randint(0,pow(2,2*nfnum)-1)
		tmp = bin(v)[2:]                     #get binary format
		inv = '0'*(2*nfnum-len(tmp))+tmp  #make tmp complete with bit number 
                #print "foundergt"+repr(foundergt)
		new_alleles = fam.getGT(inv,foundergt) #assign fam member genotypes based on inv
                #if there is missing parent, re-impute parental GT based on
                #newly generated offspring GT
                #print "new alleles"+repr(new_alleles)
                conditional_prob = {}
                if fam.missing_all != []:
                    for person in fam.missing_all:#fam.sorted_miss_persons:
                        new_alleles[person]=[0,0]
		allele_key=''
		if fam.identical_sib != []:
		    #affected siblings that are identical in terms of GT imputation
		    for affsib in fam.identical_sib:
			identical_sib_gt=[]
			for tmpsib in affsib:
			    identical_sib_gt.append(new_alleles[tmpsib])
			if identical_sib_gt == sorted(identical_sib_gt):
			    continue
			identical_sib_gt=sorted(identical_sib_gt)
			for tsib_idx,tmpsib in enumerate(affsib):
			    new_alleles[tmpsib]=identical_sib_gt[tsib_idx]
		for person in sorted(new_alleles.keys()):
		    allele_key+=''.join(map(str,new_alleles[person]))
                if allele_key not in all_ibd.keys() or sall_flag and all_ibd[allele_key][1]==0:
                    if fam.missing_all != []:
                        #print "GT_infer start"
                        #print repr(procID)+" new_alleles:"+repr(new_alleles)
                        fam.GT_infer(new_alleles,conditional_prob,infer_flag=infer_flag)
		    	#print "GT_infer end"
                    	#print "new conditional_prob"+repr(conditional_prob)
                    #number of alleles IBD among affecteds under inv
                    ibd_result=fam.cal_ibd(new_alleles,conditional_prob,pinv_pair_dict,pinv_sall,pinv_key_dict,sall_flag)
                    #print repr(procID)+" IBD finish: "+repr(ibd_result)
                    all_ibd[allele_key]=ibd_result
                    null_ibd.append(ibd_result)
                else:
                    while True:
                        try:
                            ibd_result=all_ibd[allele_key]
                            null_ibd.append(ibd_result)
                            break
                        except:
                            #in case key error, wait for the key:value to be assigned
                            pass
        except Exception as e:
            screen_output.err_out("error in null_generator {}:{}".format(procID,e))
	    raise TypeError
            break

class myProcess(Process):
    def __init__(self, procID, founderid, fam, q,all_ibd,null_ibd,\
            pinv_sall,pinv_key_dict,pinv_pair_dict,sall_flag,infer_flag=2):
        Process.__init__(self)
        self.founderid = founderid
        self.q = q
        self.fam = fam 
        self.procID = procID
	self.all_ibd = all_ibd
        self.null_ibd = null_ibd
        self.pinv_sall = pinv_sall
        self.pinv_key_dict = pinv_key_dict
        self.pinv_pair_dict = pinv_pair_dict
        self.sall_flag = sall_flag
        self.infer_flag = infer_flag
    def run(self):
        try:
            null_generator(self.procID,self.founderid, self.fam, self.q, \
                    self.all_ibd,self.null_ibd,\
                    self.pinv_sall,self.pinv_key_dict,self.pinv_pair_dict,self.sall_flag,self.infer_flag)
        except Exception as e:
            screen_output.err_out("error in: %d"%self.procID)
#############################################
class Family:
	def __init__(self,mid):
		self.fam_dict = {}              #dictionary for stating family relationship
		self.fam=[]                     #original pedigree information
		self.un=[]                      #pairs of uncle-nephew,cousin,sibpair
		self.cousin=[]
		self.sib=[]
		self.hsib=[]
                self.gp=[]
                self.offspring_pairs=[]         #one is the other's offspring
                self.affected_pairs=[]         #affected relative pairs (including offspring)
                self.ibd_total=0
                self.ibd_sall=0
                self.affected=[]
                self.nonfounder=[]
                self.parents={}
                self.mates={}
                self.offspring={}
                self.sex={}
                self.founder=[]
                self.affnf=[]                   #affected nonfounders
                self.invnf=[]                   #nonfounders that appear in inheritance vectors
                self.invnf_sorted=[]
                self.invfounder=[]
                self.fallele=[]
                self.all_ibd={}
                self.null_ibd=[]                #expected IBD (pair,all) under null hypothesis
                self.dist_s=None
                self.null_ibd_len=None
                self.pinv_sall=[]
                self.pinv_key_dict={}
                self.pinv_pair_dict={}
                self.null_mean = None
                self.null_std = None
                self.sall_null_mean = None
                self.sall_null_std = None
                #family structures; NEED further development to consider different affecteds
                self.parentsid = []
                self.gpid=None
                self.gmid=None
                self.famstruct = {}
                self.detailstruct = {}
                self.identical_sib = []
                self.info = True
                self.wt_fam = False
                self.simple = False
                self.mid = mid
                self.missing_all = []
                self.sorted_miss_persons=[]
                self.mp_freq = []
                self.conditional_prob = {}
                self.err = False
                self.rvibd = False
	def setdict(self,family):
		self.fam=family
		fam_num=len(self.fam)
		for i in range(fam_num):        #for each family member
			iid = self.fam[i][0]
                        tmp = {'parents':[],'offspring':[],'mate':[],'sex':[],'trait':[],'gt':[]}
			tmp['parents'] = self.fam[i][1:3]
			fid,mid = tmp['parents'][:]       #parents id
                        tmp['sex']=self.fam[i][3]
                        tmp['trait']= self.fam[i][4]
                        tmp['gt']=self.fam[i][5:]
			if fid != '0':
                            for pid in [m for m in [fid,mid] if m not in self.fam_dict]:
                                self.fam_dict[pid]={'parents':[],'offspring':[],'mate':[]}
                            #fam_dict={id:{parents:[],offspring:[],mate:[]},id:{},{}}
                            self.fam_dict[fid]['offspring'].append(iid)    
                            self.fam_dict[mid]['offspring'].append(iid)
                            if mid not in self.fam_dict[fid]['mate']:
                                self.fam_dict[fid]['mate'].append(mid)
                            if fid not in self.fam_dict[mid]['mate']:
                                self.fam_dict[mid]['mate'].append(fid)
                        if iid in self.fam_dict:
                            self.fam_dict[iid]['parents']=tmp['parents']
                            self.fam_dict[iid]['sex']=tmp['sex']
                            self.fam_dict[iid]['trait']=tmp['trait']
                            self.fam_dict[iid]['gt']=tmp['gt']
                        else:
			    self.fam_dict[iid]=tmp
		#print self.fam_dict
		self.set_famstruct()
        def core_struct(self,nf):
                #get fam affected structure
                total_branch=['%d'%(nf in self.affected)]
                founder_gt=[]
                self.invnf.append(nf)
                if self.fam_dict[nf]['offspring']:
                    mates=self.fam_dict[nf]['mate']
                    if len(mates)==1:
                        founder_gt.append(sorted(self.fam_dict[mates[0]]['gt'][2*self.mid:2*self.mid+2]))
                        for tmp_off in sorted(self.fam_dict[nf]['offspring'],key=lambda x: x in self.affected, reverse=True):
                            branch,f_gt=self.core_struct(tmp_off)
                            total_branch.append(branch)
                            founder_gt.append(f_gt)
                    else:
                        for mate in mates:
                            mate_branch=[]
                            mate_fgt=[sorted(self.fam_dict[mate]['gt'][2*self.mid:2*self.mid+2])]
                            for tmp_off in sorted(self.fam_dict[mate]['offspring'],key=lambda x: x in self.affected, reverse=True):
                                branch,f_gt=self.core_struct(tmp_off)
                                mate_branch.append(branch)
                                mate_fgt.append(f_gt)
                            total_branch.append(mate_branch)
                            founder_gt.append(mate_fgt)
                return total_branch,founder_gt

	def set_famstruct(self):
                # get the affecteds and nonfounder   
                self.famstruct,self.detailstruct,self.parents,self.mates,self.offspring,self.sex={},{},{},{},{},{}
                self.invnf,self.parentsid,self.affnf,self.affected,self.nonfounder,self.founder,self.missing_all,self.fallele=[],[],[],[],[],[],[],[]
                self.invfounder,self.identical_sib=[],[]
                for iid in self.fam_dict.keys():                      #for each fam member
                    if self.fam_dict[iid]['trait']==2:
                            #pick out all affected individuals
                            self.affected.append(iid)
                    if self.fam_dict[iid]['parents']!= ['0','0']:
                            self.nonfounder.append(iid)
                    elif self.fam_dict[iid]['parents']== ['0','0']:
                            self.founder.append(iid)
                    self.parents[iid]=self.fam_dict[iid]['parents']
                    self.mates[iid]=self.fam_dict[iid]['mate']
                    self.offspring[iid]=self.fam_dict[iid]['offspring']
                    self.sex[iid]=self.fam_dict[iid]['sex']
                self.missing_all=self.findall_miss()
                #remove missing members that cannot be inferred
                for tmp_p in self.fam_dict.keys():
                    if tmp_p in self.missing_all and self.fam_dict[tmp_p]['offspring']==[]:
                        self.remove(tmp_p)
                foundergt = []
                for f in self.founder:
                    foundergt += self.fam_dict[f]['gt'][2*self.mid:2*self.mid+2]
                allgt = []
                for f in self.fam_dict.keys():
                    allgt += self.fam_dict[f]['gt'][2*self.mid:2*self.mid+2]
                if not self.missing_all:
                    info_flag=False
                    for nf in self.nonfounder:
                        if not info_flag:
                            for parent in self.fam_dict[nf]['parents']:
                                tmp_gt=self.fam_dict[parent]['gt'][2*self.mid:2*self.mid+2]
                                if tmp_gt[0]!=tmp_gt[1]:
                                    info_flag=True
                                    break
                if len(set(allgt)) == 1 or not self.missing_all and not info_flag:
                    self.info = False
                if set(allgt)==set([0,1]) and self.rvibd:
                    self.wt_fam = True
                self.affnf = sorted([nf for nf in self.nonfounder if nf in self.affected])
                for affnf in self.affnf:
                    if self.fam_dict[affnf]['offspring']==[]:
                        tmpf,tmpm=self.fam_dict[affnf]['parents']
                        if len(self.fam_dict[tmpf]['mate'])==1 and len(self.fam_dict[tmpm]['mate'])==1:
                            afsib=[off for off in self.fam_dict[tmpf]['offspring'] if self.fam_dict[off]['trait']==2 and self.fam_dict[off]['offspring']==[]]
                            if len(afsib)>1 and tuple(sorted(afsib)) not in self.identical_sib:
                                self.identical_sib.append(tuple(sorted(afsib)))
                self.parentsid = [p for nf in self.nonfounder for p in self.fam_dict[nf]['parents']]
                self.parentsid = sorted(self.parentsid,key=lambda x: self.fam_dict[x]['sex'])
                first_gen=[iid for iid in self.founder if self.generation(iid)==0]
                if min([self.generation(nf) for nf in self.nonfounder]) > -3 and len(first_gen)==2:
                    self.simple=True
                    self.gpid,self.gmid=sorted(first_gen, key=lambda x: self.fam_dict[x]['sex'])
                self.famstruct['aff']=[]
                self.famstruct['info']=[]
                if len(first_gen)==2:
                    first_gen_aff=len(set(first_gen)&set(self.affected))
                    self.famstruct['aff'].append(first_gen_aff)
                    tmp_fgt_info=[]
                    for tmp_fg in sorted(first_gen,key=lambda x: x in self.affected):
                        tmp_fgt_info.append(sorted(self.fam_dict[tmp_fg]['gt'][2*self.mid:2*self.mid+2]))
                    if first_gen_aff:
                        self.famstruct['info']=tmp_fgt_info
                    else:
                        self.famstruct['info']=sorted(tmp_fgt_info)
                    for nf in sorted([iid for iid in self.fam_dict[first_gen[0]]['offspring']], key=lambda x: x in self.affected, reverse=True):
                        aff_branch, founder_info = self.core_struct(nf)
                        self.famstruct['aff'].append(aff_branch)
                        self.famstruct['info'].append(founder_info)
                else:
                    #more than one marriage
                    first_gen_parents=[]
                    for fg_p in first_gen:
                        for mate in self.fam_dict[fg_p]['mate']:
                            first_gen_p=sorted([fg_p,mate])
                            if first_gen_p not in first_gen_parents:
                                first_gen_parents.append(first_gen_p)
                    shared_parent=list(set(first_gen_parents[0])&set(first_gen_parents[1]))[0]
                    self.famstruct['info']=[sorted(self.fam_dict[shared_parent]['gt'][2*self.mid:2*self.mid+2])]
                    for first_gen_p in first_gen_parents:
                        other_parent=first_gen_p[0] if first_gen_p[1]==shared_parent else first_gen_p[1]
                        tmp_founder_info=[sorted(self.fam_dict[other_parent]['gt'][2*self.mid:2*self.mid+2])]
                        branch=[('%d'%(shared_parent in self.affected),'%d'%(other_parent in self.affected))]
                        for nf in sorted([iid for iid in self.fam_dict[first_gen_p[0]]['offspring'] if sorted(self.fam_dict[iid]['parents'])==first_gen_p],key=lambda x: x in self.affected, reverse=True):
                            aff_branch, founder_info = self.core_struct(nf)
                            branch.append(aff_branch)
                            tmp_founder_info.append(founder_info)
                        self.famstruct['aff'].append(branch)
                        self.famstruct['info'].append(tmp_founder_info)
                self.invnf_sorted = sorted(self.invnf, key=lambda x: self.generation(x), reverse=True)
                for inf in self.invnf:
                    for parent in self.fam_dict[inf]['parents']:
                        if parent in self.founder and parent not in self.invfounder:
                            self.invfounder.append(parent)
                self.invfounder=sorted(self.invfounder,key=lambda x:self.generation(x),reverse=True)
                for f in self.invfounder:
                    self.fallele+=['{0}0'.format(f),'{0}1'.format(f)]
 
        def clean(self):
            self.offspring_pairs, self.affected_pairs = [], []
            self.null_ibd,self.pinv_sall,self.fallele = [], [], []
            self.pinv_key_dict,self.pinv_pair_dict,self.all_ibd={},{},{}

	def firstparent(self,parent,gp,gm):         #make the first parent the important one 
		#if mother is the offspring of grandparents
                if self.fam_dict[parent[1]]['parents'] == [gp,gm]:        
			parent.reverse()            #put the important parent in the first

	def dcousin(self,ap,bp,gp,gm):                          #check if they are cousin
                if ap != ['0','0'] and ap != [gp,gm] and bp != ['0','0'] and bp !=[gp,gm] \
				and set(ap)&set(bp)==set([]):               #if the affecteds are 1st cousin
			return True
		else:
			return False

	def dun(self,aid,bid,ap,bp,gp,gm):              #check if it is uncle-nephew pair
		i = False
                if ap !=['0','0'] and bp !=['0','0'] and bid not in ap and aid not in bp:
			if ap == [gp,gm] and bp != [gp,gm] or ap != [gp,gm] and bp == [gp,gm]:
				i = True
		return i		

	def dsib(self,ap,bp):                       		#check if it is sibpair
		i = False
		if ap == bp and ap != ['0','0']:
			i = True
		return i

	def dhsib(self,ap,bp):
		i = False
		if ap != bp and len(set(ap)&set(bp)) > 0:
			i = True
		return i

        def dgrandpa(self,aid,bid,ap,bp):
                if self.generation(aid)==self.generation(bid)+2:
                        if set(bp)&set(self.fam_dict[aid]['offspring'])!=set([]):
                                return 1
                elif self.generation(aid)==self.generation(bid)-2:
                        if set(ap)&set(self.fam_dict[bid]['offspring'])!=set([]):
                                return 2
                else:
                        return 0

        def dchildren(self,aid,bid):
                if aid in self.fam_dict[bid]['parents']:
                        return 1
                elif bid in self.fam_dict[aid]['parents']:
                        return 2
                else:
                        return 0

        def doffspring(self,aid,bid):
                gen_a=self.generation(aid)
                gen_b=self.generation(bid)
                correct=False
                if gen_a==gen_b:
                        return False
                else:
                        diff_gen=abs(gen_a-gen_b)
                        older=aid if gen_a>gen_b else bid
                        younger=aid if older==bid else bid
                        if younger in self.fam_dict[older]['offspring']:
                                correct=True
                        if not correct:
                                for off in self.fam_dict[older]['offspring']:
                                        if self.doffspring(off,younger):
                                                correct=True
                                                break
                        return correct

        def generation(self,iid):
                gen=0
                local_id=iid
                while(self.fam_dict[local_id]['parents']!=['0','0']):
                    gen-=1
                    nf_parent=[parent for parent in self.fam_dict[local_id]['parents'] if self.fam_dict[parent]['parents']!=['0','0']]
                    if len(nf_parent)==0:
                        break
                    elif len(nf_parent)==1:
                        local_id=nf_parent[0]
                    else:
                        raise TypeError
                        break
                if gen==0:
                    mate=self.fam_dict[local_id]['mate'][0]
                    if self.fam_dict[mate]['parents']!=['0','0']:
                        gen=self.generation(mate)
                return gen
		    
        
#classify affecteds and calculate the IBD allele number
        def classify_affect(self):                   #mid: marker index
                if self.simple:
                    affected_num = len(self.affected)
                    for i in range(affected_num):
                            for j in range(i+1,affected_num):   #for each pair of affecteds
                                    aid = self.affected[i]
                                    bid = self.affected[j]
                                    ap = self.fam_dict[aid]['parents']
                                    bp = self.fam_dict[bid]['parents']
                                    if self.dcousin(ap,bp,self.gpid,self.gmid):
                                        # if a and b are cousin
                                            self.cousin.append([aid,bid])
                                    elif self.dun(aid,bid,ap,bp,self.gpid,self.gmid):
                                        # uncle-nephew
                                            if bp == [self.gpid,self.gmid]:   #put uncle at first
                                                    self.un.append([bid,aid])
                                            else:
                                                    self.un.append([aid,bid])
                                    elif self.dsib(ap,bp):                              # sibpair
                                            self.sib.append([aid,bid])
                                    elif self.dhsib(ap,bp):
                                            self.hsib.append([aid,bid])
                                    tmp_gp=self.dgrandpa(aid,bid,ap,bp)
                                    if tmp_gp==1:
                                            self.gp.append([aid,bid])
                                    elif tmp_gp==2:
                                            self.gp.append([bid,aid])
                miss_persons=[]
                for member in self.missing_all:
                    offmiss_count=0
                    for tmp_off in self.fam_dict[member]['offspring']:
                        if tmp_off in self.missing_all:
                            offmiss_count+=1
                    miss_persons.append((member,offmiss_count))
                        #sort miss_persons so that person with sequenced offspring can\
                        #be inferred first
                self.sorted_miss_persons=sorted(miss_persons,key=lambda x: (self.generation(x[0]),x[1]))
                if self.sorted_miss_persons != []:
                    self.GT_infer()
                #print self.conditional_prob
                #remove useless individuals
                rm_flag=False
                rm_genotyped_person=False
                if not self.err and self.sorted_miss_persons != []:
                    for m_person in self.sorted_miss_persons:
                        if m_person[0] not in self.conditional_prob['~combined'][1]:
                            self.remove(m_person[0])
                            rm_flag=True
                    try:
                        if self.conditional_prob['~combined'][0]=={}:
                            screen_output.err_out('conditional_prob is NULL')
                            self.err=True
                    except KeyError:
                        screen_output.err_out('conditional_prob no key')
                        self.err = True
                if rm_flag:
                    self.set_famstruct()
                self.offspring_pairs=[]
                for idx,aid in enumerate(self.affected):
                    for bid in self.affected[idx+1:]:
                        if aid in self.founder and bid in self.founder or aid in self.fam_dict[bid]['mate']:
                            continue
                        dchildren=self.dchildren(aid,bid)
                        if dchildren:
                            if dchildren==1:
                                self.offspring_pairs.append([bid,aid])
                            else:
                                self.offspring_pairs.append([aid,bid])
                            self.affected_pairs.append([aid,bid])
                        else:
                            if aid in self.founder or bid in self.founder:
                                if not self.doffspring(aid,bid):
                                    continue
                            self.affected_pairs.append([aid,bid])
                if len(self.affected_pairs)==0:
                    screen_output.err_out('no affected relatives pairs')
                    self.err = True       
	
	def remove(self,iid):
	        #remove useless individuals
	    if self.fam_dict.has_key(iid):
		#print "removing {}".format(iid)
		if iid in self.missing_all:
		    for m_person in self.sorted_miss_persons:
			if iid in m_person:
		    	    self.sorted_miss_persons.remove(m_person)
		    self.missing_all.remove(iid)
		for tmp_p in self.fam_dict.keys():
		    if iid in self.fam_dict[tmp_p]['mate']:
			self.fam_dict[tmp_p]['mate'].remove(iid)
		    if iid in self.fam_dict[tmp_p]['offspring']:
			self.fam_dict[tmp_p]['offspring'].remove(iid)
		    if iid in self.fam_dict[tmp_p]['parents']:
			idx=self.fam_dict[tmp_p]['parents'].index(iid)
			self.fam_dict[tmp_p]['parents'][idx]='0'
			
		try:
		    self.affected.remove(iid)
		except:
		    pass
		try:
		    self.founder.remove(iid)
		except:
		    pass
		try:
		    self.nonfounder.remove(iid)
		except:
		    pass
		try:
		    self.affnf.remove(iid)
		except:
		    pass
		try:
		    self.invnf.remove(iid)
		except:
		    pass
		try:
		    self.parentsid.remove(iid)
		except:
		    pass
                for afidx,affsibs in enumerate(self.identical_sib):
                    if iid in affsibs:
                        list_affsibs=list(affsibs)
                        list_affsibs.remove(iid)
                        if len(list_affsibs)>1:
                            self.identical_sib[afidx]=tuple(list_affsibs)
                        else:
                            self.identical_sib.remove(affsibs)
		#for off in self.fam_dict[iid]['offspring']:
		#    self.remove(off)
		try:
		    del self.fam_dict[iid]
		except:
		    pass
                self.missing_all=self.findall_miss()
		#remove missing members that cannot be inferred
		for tmp_p in self.fam_dict.keys():
		    if tmp_p in self.missing_all and self.fam_dict[tmp_p]['offspring']==[]:
			try:
			    self.remove(tmp_p)
			except:
			    pass

        def findall_miss(self):
            #find all individuals that have missing GT
            all_gt=[]
            all_members=[]
            for iid in self.fam_dict.keys():
                all_members.append(iid)
                all_gt+=self.fam_dict[iid]['gt'][2*self.mid:2*self.mid+2]
            if 0 in all_gt and len(set(all_gt))>1:
                #there are missing individuals
                miss=[all_members[int(idx/2)] for idx,a in enumerate(all_gt) if a==0]
                return list(set(miss))
            return []

        def conv2(self,dic):
            tmp_dic=copy.copy(dic)
            key_list=sorted(tmp_dic.keys(),reverse=True)
            first_key=key_list.pop(0)
            first_list=tmp_dic.pop(first_key,None)
            return_gt=[]
            if tmp_dic.keys()==[]:
                #last one
                for tmp_ele in first_list:
                        return_gt.append([tmp_ele])
            else:
                pre_gt=self.conv2(tmp_dic)
                for gt in pre_gt:
                    for tmp_ele in first_list:
                        combined_gt=gt+[tmp_ele]
                        return_gt.append(combined_gt)
            return return_gt

        def infer_engine(self,couple,alleles=[],possible_GT=[],hsib_flag=False):
            #genotype inference
	    num_to_letter = dict(enumerate(string.ascii_lowercase, 1))
            father = couple[0] if self.fam_dict[couple[0]]['sex']==1 else couple[1]
            couple.remove(father)
            mother = couple[0]
            pgt=alleles[father]+alleles[mother]
            tmp_offgt=[]
            mp_offgt=[]
            mp_pgt=[]
            if 0 in pgt[:2]:
                mp_pgt.append("None")
            else:
		mp_pgt.append(num_to_letter[int(pgt[0])]+num_to_letter[int(pgt[1])])
            if 0 in pgt[2:]:
                mp_pgt.append("None")
            else:
		mp_pgt.append(num_to_letter[int(pgt[2])]+num_to_letter[int(pgt[3])])
            if self.mp_freq == []:
                raise Exception('marker frequency info missing')
            miss_off_gt={}
            if len(self.fam_dict[father]['mate'])==1:
                offspring=sorted(self.fam_dict[father]['offspring'])
            else:
                offspring=sorted(self.fam_dict[mother]['offspring'])
            offspring_idx=[]
            combined_id=[]
            for off in offspring:
                try:
                    cond_P=possible_GT[off][0]
                    tmp_ids=possible_GT[off][1]
                    offspring_idx.append(tmp_ids.index(off))
                    miss_off_gt[off]=cond_P
                    combined_id+=tmp_ids
                except:
                    offspring_idx.append(0)
                    tmp_gt=alleles[off]
                    if tmp_gt==[None,None]:
                        #if this offspring is missing in GT, and doesn't have imputation, then ignore
                        continue
                    key=num_to_letter[int(tmp_gt[0])]+num_to_letter[int(tmp_gt[1])]
                    combined_id+=[off]
                    miss_off_gt[off]=[key]
            #print "miss_off_gt:"+repr(miss_off_gt)
            return_gt =self.conv(miss_off_gt)
            combined_possible_GT=[]
            for rgt in return_gt:
                mp_offgt=[]
                for off_idx,t_rgt in enumerate(rgt):
                    id_idx=offspring_idx[off_idx]
                    mp_offgt.append(t_rgt[2*id_idx:2*id_idx+2])
                tmp_possible_GT=cmissingparents.apply(self.mp_freq,mp_offgt,mp_pgt)
                combined_offgt=''
                for offgt in rgt:
                    combined_offgt+=offgt
                for pgt in tmp_possible_GT:
                    combined_gt=combined_offgt+pgt
                    combined_possible_GT.append(combined_gt)
            return [combined_possible_GT,combined_id+[father,mother]]
        
        def hsib_infer(self,dic,shared,infer_persons):
            shared_pos=0 if self.fam_dict[shared]['sex']==1 else 1
            mates=sorted(dic.keys())
            pre_combined_gt=self.conv2(dic)
	    #filter out those cases that have conflict genotypes on shared partner
            for idx,tmp_gt in enumerate(pre_combined_gt):
                shared_gt=[]
                for pair_gt in tmp_gt:
		    if len(pair_gt)>4:
			pair_gt=pair_gt[-4:]
                    shared_gt.append(tuple(pair_gt[2*shared_pos:2*shared_pos+2]))
                if len(set(shared_gt))>1:
                    pre_combined_gt[idx]=None
            new_combined_gt=[x for x in pre_combined_gt if x != None]
            #print "new_combined_gt:"+repr(new_combined_gt)
            possible_GT=[]
            for gt in new_combined_gt:
                    tmp_pgt=gt[0][-4:]
                    key=tmp_pgt[2*shared_pos:2*shared_pos+2]
                    for tmp_gt in gt:
                        parent_gt=tmp_gt[-4:]
                        key+=parent_gt[2*(1-shared_pos):2*(1-shared_pos)+2]
                        if len(tmp_gt)>4:
                        #if there are more people
                            key+=tmp_gt[:-4]
                    possible_GT.append(key)
            combined_id=[shared]
            for mate in mates:
                combined_id.append(mate)
                if len(infer_persons[mate])>2:
                    combined_id+=infer_persons[mate][:-2]
            return [possible_GT,combined_id]

        def GT_infer(self,alleles=[],conditional_prob=[],infer_flag=2):
            #Dealing with missing relatives of sequenced individuals
	    initial_flag=False
            if conditional_prob==[]:
                conditional_prob=self.conditional_prob
            if alleles == []:
		initial_flag=True
                alleles={}
                for iid in self.fam_dict.keys():
                    alleles[iid] = self.fam_dict[iid]['gt'][2*self.mid:2*self.mid+2]
            sorted_miss_persons=[]
            for tmp_miss in self.sorted_miss_persons:
                sorted_miss_persons.append(tmp_miss[0])
            tmp_dic,miss_gt_id=cmissing_infer.apply(self.mp_freq,self.invnf,self.founder,self.parents,self.mates,self.offspring,self.sex,sorted_miss_persons,alleles,infer_flag)
            cut_size=50
            cut_prob=1
            if len(tmp_dic.keys())>cut_size:
                key_gt, cut_prob=sorted([(k,v) for k,v in tmp_dic.iteritems()],key=lambda x: x[1], reverse=True)[cut_size]
            if cut_prob==1:
                cut_prob=1E-6
            for key_gt in tmp_dic.keys():
                #limit posterior imputation prob by removing unlikely imputed genotypes
                if tmp_dic[key_gt]<cut_prob:
                    del tmp_dic[key_gt]
            sum_prob=sum([tmp_dic[key] for key in tmp_dic.keys()])
            if sum_prob != 1.0:
                for tmp_key in tmp_dic.keys():
                    tmp_dic[tmp_key]=tmp_dic[tmp_key]/sum_prob
            #print "GT_infer finished"
            conditional_prob['~combined']=[tmp_dic,miss_gt_id]

        def conv(self,dic):
            #combine the inheritance vectors of nonfounders
            d = copy.copy(dic)
            k=sorted(d.keys())
            result = []
            if len(k) > 2:
                first_key=k[0]
                pred=d[first_key]
                del d[first_key]
                con = self.conv(d)
            elif len(k) == 2:
                i,j = tuple(k)
                pred = d[i]
                con = d[j]
            elif len(k) == 1:
                return [[ele] for ele in d[k[0]]]
            for ele1 in pred:
                for ele2 in con:
                    tmp_result=[ele1]
                    if isinstance(ele2,list):
                        tmp_result.extend(ele2)
                    else:
                        tmp_result.append(ele2)
                    result.append(tmp_result)
            return result

        def assign_allele(self,v):
            fam_allele={}
            ##assign alleles from nonfounders as founder alleles
            for mid in self.invnf_sorted:
                idx=self.invnf.index(mid)
                fam_allele[mid]=[None,None]
                mv=[int(x) for x in v[2*idx:2*idx+2]]
                mp = self.fam_dict[mid]['parents']
                try:
                    fam_allele[mid][0]=fam_allele[mp[0]][mv[0]]
                except:
                    fam_allele[mid][0]='{0}{1}'.format(mp[0],mv[0])
                try:
                    fam_allele[mid][1]=fam_allele[mp[1]][mv[1]]
                except:
                    fam_allele[mid][1]='{0}{1}'.format(mp[1],mv[1])
            for founder in self.founder:
                fam_allele[founder]=['{0}0'.format(founder),'{0}1'.format(founder)]
            return fam_allele

        def sall(self,v,fam_allele,fmarker):
            ##calculate S-all based on given inheritance vector#####
            h = []           #collections of alleles from affected individuals
            try:
                self.affnf[0]+1
                s_all = sall_cpp.apply(self.affected,fam_allele,self.fallele,fmarker,self.rvibd)
            except TypeError:
                #code affected if not int
                code_dic={}
                code_int=0
                for tmp_af in self.affected:
                    code_int+=1
                    code_dic[tmp_af]=code_int
                tmp_aff=[code_dic[tmp_af] for tmp_af in self.affected]
                tmp_fam_allele={}
                for tmp_key in self.affected:
                    tmp_fam_allele[code_dic[tmp_key]]=fam_allele[tmp_key]
                if fmarker==[]:
                    s_all = sall_cpp.apply(tmp_aff,tmp_fam_allele,self.fallele,fmarker,False)
                else:
                    s_all = sall_cpp.apply(tmp_aff,tmp_fam_allele,self.fallele,fmarker,self.rvibd)
            return s_all
	
        def getGT(self,v,foundergt):
            #assign nonfounders' genotypes based on inv and founder genotypes
                alleles=copy.copy(foundergt)
                for mid in sorted(self.nonfounder, key=lambda x: self.generation(x), reverse=True):
                        alleles[mid]=[]
                        parents=self.fam_dict[mid]['parents']
                        v_id=self.nonfounder.index(mid)
                        mv=v[2*v_id:2*v_id+2]
                        alleles[mid].append(alleles[parents[0]][int(mv[0])])  #parental
                        alleles[mid].append(alleles[parents[1]][int(mv[1])])  #maternal
                for mid in alleles.keys():
                        alleles[mid]=sorted(alleles[mid])
                return alleles

        def cal_sall(self,pinv,fmarker,pinv_sall=[],pinv_key_dict={},tmp_alleles={},t_alleles={},rvibd=False):
            tmp_sall=0
            ###NPL-all
            if rvibd: #if RV only, sall is dependent on founder alleles and post inv
                aff_rv_alleles=[]  #inheritance code for RVs in affecteds
                dup_person={}
                multi_a={}
                for af in self.affected:
                    if tmp_alleles[af] != [1,1]:
                        for allele_inv in t_alleles[af]:
                            if allele_inv in aff_rv_alleles:
                                aff_rv_alleles.append(allele_inv)
                                continue
                            founder_id = allele_inv[:-1]
                            tmp_inv = int(allele_inv[-1:])
                            marker = tmp_alleles[founder_id][tmp_inv]
                            if marker != 1:
                                aff_rv_alleles.append(allele_inv)
                tmp_code=[]
                for founder_a in set(aff_rv_alleles):
                        count=aff_rv_alleles.count(founder_a)
                        if count>1:
                            tmp_code.append(count)
                            multi_a[founder_a]=[count,tmp_code.count(count)]
                if tmp_code==[]:
                    tmp_sall=1
                    return tmp_sall
                if len(tmp_code)==1:
                    pinv_key=tmp_code[0]
                else:
                    tmp_code=sorted(tmp_code)
                    multi_a_keys=multi_a.keys()
                    for afid in self.affected:
                        if t_alleles[afid][0] in multi_a_keys and t_alleles[afid][1] in multi_a_keys:
                        #if len(set(t_alleles[afid])&set(multi_a.keys()))==2:
                            dp_key=tuple(sorted(['-'.join(map(str,multi_a[t_alleles[afid][0]])),'-'.join(map(str,multi_a[t_alleles[afid][1]]))]))
                            try:
                                dup_person[dp_key]+=1
                            except:
                                dup_person[dp_key]=1
                    pinv_key=tuple([tuple(tmp_code),tuple(sorted(dup_person.iteritems()))])
                try:
                    tmp_sall = pinv_key_dict[pinv_key]
                except:
                    #print self.fallele, t_alleles, self.affected
                    #print tmp_alleles
                    #print "sall"
                    tmp_sall=self.sall(pinv,t_alleles,fmarker)
                    pinv_key_dict[pinv_key]=tmp_sall
                #print "#",tmp_sall,pinv_key
                #print "#"*10
            else: #if including all variants, sall is dependent only on post inv
                int_pinv=int(pinv,2)
                if pinv_sall[int_pinv] != 0.0:
                    tmp_sall=pinv_sall[int_pinv]
                else:
                    aff_allele_pool=[]
                    tmp_code=[]
                    dup_person={}
                    multi_a={}
                    for afid in self.affected:
                            aff_allele_pool.extend(t_alleles[afid])
                    for founder_a in list(set(aff_allele_pool)):
                        fa_count=aff_allele_pool.count(founder_a)
                        if fa_count>1:
                            tmp_code.append(fa_count)
                            multi_a[founder_a]=[fa_count,tmp_code.count(fa_count)]
                    for afid in self.affected:
                        if len(set(t_alleles[afid])&set(multi_a.keys()))==2:
                            dp_key=tuple(sorted(['-'.join(map(str,multi_a[t_alleles[afid][0]])),'-'.join(map(str,multi_a[t_alleles[afid][1]]))]))
                            try:
                                dup_person[dp_key]+=1
                            except:
                                dup_person[dp_key]=1
                    pinv_key=tuple([tuple(sorted(tmp_code)),tuple(sorted(dup_person.iteritems()))])
                    if pinv_key in pinv_key_dict:
                        pinv_sall[int_pinv]=pinv_key_dict[pinv_key]
                        tmp_sall=pinv_key_dict[pinv_key]
                    else:
                        tmp_sall=self.sall(pinv,t_alleles,fmarker)
                        pinv_key_dict[pinv_key]=tmp_sall
                        pinv_sall[int_pinv]=tmp_sall
            return tmp_sall

        def cal_ibd(self,alleles={},conditional_prob=[],pinv_pair_dict={},pinv_sall=[],pinv_key_dict={},sall_flag=False):
                #calculate both NPL-pairs and NPL-all
                tmp_alleles=copy.copy(alleles)
                if conditional_prob==[]:
                    conditional_prob=self.conditional_prob
                if tmp_alleles=={}:
                    for iid in self.fam_dict.keys():
                        tmp_alleles[iid]=self.fam_dict[iid]['gt'][2*self.mid:2*self.mid+2]
                if pinv_pair_dict == {}:
                    pinv_pair_dict = self.pinv_pair_dict
                if pinv_sall == []:
                    pinv_sall = self.pinv_sall
                    if self.pinv_sall == []:
                        self.pinv_sall = Array(c_double,[0 for x in xrange(pow(2,2*len(self.invnf)))])
                        pinv_sall=self.pinv_sall
                if pinv_key_dict == {}:
                    pinv_key_dict = self.pinv_key_dict
                tmp_dic={}
                gt_id=[]
                output_ibd=0
                output_sall=0
                if self.sorted_miss_persons != []:
                    tmp_dic=conditional_prob['~combined'][0]
                    gt_id=conditional_prob['~combined'][1]
                else:
                    tmp_dic={'0000':1}  #dummy dic when no missing
                #print "prob_gt:%d"%(len(tmp_dic.keys()))
                #t0=time.time()
                for tmp_gt,prob in tmp_dic.iteritems():
                    total_ibd = 0
                    total_sall = 0
                    #print tmp_gt,prob
                    #fill out missing genotypes based on inference
                    for m_person in self.sorted_miss_persons:
                        idx=gt_id.index(m_person[0])
                        tmp_alleles[m_person[0]]=[ord(x)-96 for x in tmp_gt[2*idx:2*idx+2]]
                    if self.rvibd:
                        cal_flag=False
                        affmarker=[]
                        if sall_flag:
                            for afid in self.affected:
                                affmarker+=tmp_alleles[afid]
                            for uniq_marker in list(set(affmarker)):
                                if uniq_marker!=1 and affmarker.count(uniq_marker)>1:
                                    cal_flag=True
                                    break
                        else:
                            aff_ind=list(set([tmp_aff for pair in self.affected_pairs for tmp_aff in pair]))
                            for afid in aff_ind:
                                affmarker+=tmp_alleles[afid]
                            for uniq_marker in list(set(affmarker)):
                                if uniq_marker!=1 and affmarker.count(uniq_marker)>1:
                                    tmp_share_pattern={}
                                    for afid in aff_ind:
                                        tmp_share_pattern[afid]=[True,True] if uniq_marker in tmp_alleles[afid] else [False,False]
                                        if self.generation(afid)<-1 and afid not in self.founder and uniq_marker in tmp_alleles[afid]:
                                            nf_parent=[tmp_parent for tmp_parent in self.parents[afid] if tmp_parent not in self.founder][0]
                                            f_parent=self.parents[afid][1-self.parents[afid].index(nf_parent)]
                                            if uniq_marker not in tmp_alleles[nf_parent] and f_parent not in self.affected:
                                                tmp_share_pattern[afid][1]=False
                                    for tmp_pair in self.affected_pairs:
                                        if tmp_share_pattern[tmp_pair[0]][0] and tmp_share_pattern[tmp_pair[1]][0]:
                                            if self.parents[tmp_pair[0]]==self.parents[tmp_pair[1]]:
                                                cal_flag=True
                                                break
                                            elif tmp_share_pattern[tmp_pair[0]][1] or tmp_share_pattern[tmp_pair[1]][1]:
                                                cal_flag=True
                                                break
                                    if cal_flag:
                                        break
                        #print "cal_flag:"+repr(cal_flag)
                        if not cal_flag:
                            total_ibd = 0
                            if sall_flag:
                                total_sall = 1
                                output_sall+=total_sall*prob
                            continue
                    fmarker=[]
                    for f in self.invfounder:
                        fmarker+=tmp_alleles[f]
                    if self.simple:
                        #NPL-pair
                        #print self.sib, self.cousin, self.un
                        gpid,gmid=self.gpid,self.gmid
                        ibd_n=0        #number of alleles IBD among affected individuals        
                        for aid,bid in self.offspring_pairs:
                            if not self.rvibd:
                                total_ibd+=1
                            else:
                                if tmp_alleles[aid]==[1,1] or tmp_alleles[bid]==[1,1]:
                                    continue
                                else:
                                    other_parent=[tmp_p for tmp_p in self.fam_dict[aid]['parents'] if tmp_p != bid][0]
                                    genotype = tmp_alleles[aid]+tmp_alleles[bid]+tmp_alleles[other_parent]
                                    total_ibd+=ibd_module.child_ibd(genotype)
                        for aid,bid in self.sib:
                            ap = self.fam_dict[aid]['parents']
                            genotype = tmp_alleles[ap[0]]+tmp_alleles[ap[1]]+tmp_alleles[aid]+tmp_alleles[bid]
                            ibd = ibd_module.sib_ibd(genotype,rv_flag=self.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.cousin:
                            ap = copy.copy(self.fam_dict[aid]['parents'])
                            bp = copy.copy(self.fam_dict[bid]['parents'])
                            self.firstparent(ap,gpid,gmid)
                            self.firstparent(bp,gpid,gmid)
                            genotype = tmp_alleles[gpid]+tmp_alleles[gmid]+tmp_alleles[ap[0]]+\
                                    tmp_alleles[ap[1]]+tmp_alleles[aid]+tmp_alleles[bp[0]]+\
                                    tmp_alleles[bp[1]]+tmp_alleles[bid]
                            ibd = ibd_module.cousin_ibd(genotype,rv_flag=self.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.hsib:
                            ap = copy.copy(self.fam_dict[aid]['parents'])
                            bp = copy.copy(self.fam_dict[bid]['parents'])
                            sp = list(set(ap)&set(bp))[0]                #shared parent between half-sibs
                            ap.remove(sp)
                            bp.remove(sp)
                            genotype = tmp_alleles[sp]+tmp_alleles[ap[0]]+tmp_alleles[aid]+\
                            tmp_alleles[bp[0]]+tmp_alleles[bid]               #[shared parent, non-share parent, kid]
                            ibd=ibd_module.hsib_ibd(genotype,rv_flag=self.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.un:      #(uncle,nephew)
                            ap = copy.copy(self.fam_dict[aid]['parents'])
                            bp = copy.copy(self.fam_dict[bid]['parents'])
                            self.firstparent(bp,gpid,gmid)
                            genotype = tmp_alleles[gpid]+tmp_alleles[gmid]+tmp_alleles[aid]+\
                                    tmp_alleles[bp[0]]+tmp_alleles[bp[1]]+tmp_alleles[bid]
                            ibd = ibd_module.un_ibd(genotype,rv_flag=self.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.gp:      #(grandpa,grandson)
                            bp = copy.copy(self.fam_dict[bid]['parents'])
                            self.firstparent(bp,gpid,gmid)
                            grand_other=gmid if aid==gpid else gpid
                            genotype = tmp_alleles[aid]+tmp_alleles[grand_other]+\
                                    tmp_alleles[bp[0]]+tmp_alleles[bp[1]]+tmp_alleles[bid]
                            ibd = ibd_module.gp_ibd(genotype,rv_flag=self.rvibd)
                            total_ibd+=ibd
                        if sall_flag:
                            ###NPL-all
                            #get posterior inheritance vectors for probable allele configurations
                            postinv=cpostInv.apply(self.invnf, self.founder, self.parents, tmp_alleles, self.rvibd)
                            len_postinv=len(postinv)
                            #print fmarker
                            #print "postinv:"+repr(len_postinv)+" pinv_key_dict:"+repr(len(pinv_key_dict.keys()))
                            postinv_list=set(postinv) if self.rvibd else postinv
                            for pinv in postinv_list:
                                t_alleles=self.assign_allele(pinv)
                                t_sall=self.cal_sall(pinv,fmarker,pinv_sall,pinv_key_dict,tmp_alleles,t_alleles,rvibd=self.rvibd)
                                if self.rvibd:
                                    total_sall+=t_sall*postinv.count(pinv)
                                else:
                                    total_sall+=t_sall
                            total_sall=total_sall/len_postinv
                    else:
                        if self.rvibd:
                            founder_rv_alleles=[]
                            for founder in self.founder:
                                for tmp_i in range(2):
                                    if tmp_alleles[founder][tmp_i] not in [0,1]:
                                        #alt allele, i.e. RV
                                        founder_rv_alleles.append('{}{}'.format(founder,tmp_i))
                        #get posterior inheritance vectors for probable allele configurations
                        postinv=cpostInv.apply(self.invnf, self.founder, self.parents, tmp_alleles, self.rvibd)
                        len_postinv=len(postinv)
                        #print "postinv:"+repr(len_postinv)+" pinv_pair:"+repr(len(pinv_pair_dict.keys()))
                        postinv_list=set(postinv) if self.rvibd else postinv
                        for pinv in postinv_list:
                            t_alleles={}
                            ###NPL-pairs
                            int_pinv=int(pinv,2)
                            shared_alleles=[]
                            t_ibd=0
                            if int_pinv in pinv_pair_dict:
                                shared_alleles=pinv_pair_dict[int_pinv]
                            else:
                                t_alleles=self.assign_allele(pinv)
                                for af_pair in self.affected_pairs:
                                    for tt_a in t_alleles[af_pair[0]]:
                                        if tt_a in t_alleles[af_pair[1]]:
                                            shared_alleles.append(tt_a)
                                    #tmp_ibd=set(t_alleles[af_pair[0]])&set(t_alleles[af_pair[1]])
                                    #shared_alleles.extend(list(tmp_ibd))
                                pinv_pair_dict[int_pinv]=shared_alleles
                            if self.rvibd:
                                t_ibd=len([tmp_sh_a for tmp_sh_a in shared_alleles if tmp_sh_a in founder_rv_alleles])
                                total_ibd+=t_ibd*postinv.count(pinv)
                            else:
                                t_ibd=len(shared_alleles)
                                total_ibd+=t_ibd
                            if sall_flag:
                                ###NPL-all
                                if t_alleles=={}:
                                    t_alleles=self.assign_allele(pinv)
                                t_sall=self.cal_sall(pinv,fmarker,pinv_sall,pinv_key_dict,tmp_alleles,t_alleles,self.rvibd)
                                if self.rvibd:
                                    total_sall+=t_sall*postinv.count(pinv)
                                else:
                                    total_sall+=t_sall
                        total_ibd=total_ibd/len_postinv
                        total_sall=total_sall/len_postinv
                    output_sall+=total_sall*prob
                    output_ibd+=total_ibd*prob
                output_ibd=float('%.9f'%(output_ibd))
                output_sall=float('%.9f'%(output_sall))
                return output_ibd,output_sall

	def std(self,ibd,mean=0,sample_flag=True):
		total = 0
                ibd_size = len(ibd)
		if mean == 0:
			mean = sum(ibd)/ibd_size
		for element in ibd:
			total+=(element-mean)*(element-mean)
		if sample_flag:
		    std = math.sqrt(total/(ibd_size-1))
		if not sample_flag:
		    std = math.sqrt(total/ibd_size)
		return std

        def nullibd(self,rep,n_jobs,sall_flag=False,infer_flag=2,simple=False):
            #calculate expected mean and std for IBD under H0 
            screen_output.run_out("calculating null distribution...")
            #the number of nonfounders that should be included in inheritance vector
            founderid = self.founder
            ###parallel processing###
            inqueue = multiprocessing.Queue()
            all_ibd = manager.dict(self.all_ibd)
            null_ibd = manager.list(self.null_ibd)
            if self.pinv_sall==[]:
                pinv_sall = Array(c_double,[0 for x in xrange(pow(2,2*len(self.invnf)))])
            else:
                pinv_sall=(c_double*len(self.pinv_sall))(*self.pinv_sall)
            pinv_key_dict = manager.dict(self.pinv_key_dict)
            pinv_pair_dict = manager.dict(self.pinv_pair_dict)
            for i in xrange(rep):
                inqueue.put(i)
            #start n seperate processes
            procs = []
            for proc in range(n_jobs):
                p = myProcess(proc,founderid,self,inqueue,all_ibd,\
                        null_ibd,pinv_sall,pinv_key_dict,pinv_pair_dict,sall_flag,infer_flag)
                p.start()
                procs.append(p)
                inqueue.put(None)
            TIMEOUT = 3600 #10000 if sall_flag else 3600 
            start = time.time()
            last_flag=0
            while time.time() - start <= TIMEOUT:
                if any(p.is_alive() for p in procs):
                    if last_flag>10:
                        for p in procs:
                            p.terminate()
                            p.join()
                        break
                    elif len(null_ibd)==rep-1:
                        last_flag+=1
                    else:
                        time.sleep(.1)  # Just to avoid hogging the CPU
                else:
                    break
            else:
                try:
                    #print("timed out, killing all processes")
                    for p in procs:
                        p.terminate()
                        p.join()
                except:
                    pass
            while not inqueue.empty():
                inqueue.get()
            self.null_ibd = null_ibd
            self.pinv_pair_dict = pinv_pair_dict
            self.all_ibd = all_ibd
            if sall_flag:
                self.pinv_sall = [x for x in pinv_sall]
                self.pinv_key_dict = pinv_key_dict
            if not simple:
                #expectation and standard deviation
                ibd_pair=[x[0] for x in self.null_ibd]
                self.null_mean = float('%.9f'%(sum(ibd_pair)/len(ibd_pair)))
                self.null_std = self.std(ibd_pair,self.null_mean)
                if sall_flag:
                    ibd_all=[x[1] for x in self.null_ibd if x[1] != 0]
                    self.sall_null_mean = float('%.9f'%(sum(ibd_all)/len(ibd_all)))
                    self.sall_null_std = self.std(ibd_all,self.sall_null_mean)

	def null_perfect(self,n_jobs,perfect_max,update_null=True):
	#null distribution of IBDs under perfect data
	    screen_output.run_out("calculating null distribution using perfect data approximation...")
	    p = Pool(processes=n_jobs)
	    inv_list=[]
	    rep=pow(2,2*len(self.invnf))
	    #print rep
	    if rep>perfect_max:
		inv_list=random.sample(xrange(rep),perfect_max)
	    else:
		inv_list=range(rep)
            self.pinv_sall = []
            chunk_size=len(inv_list)//n_jobs if len(inv_list)>=n_jobs else 1
            result = p.imap(null_inv, product(inv_list,[self]), chunksize=chunk_size)
            p.close()
            p.join()
            null_ibd=[]
            invs=[]
            for r in result:
                null_ibd.append((float('%.9f'%r[0]),float('%.9f'%r[1])))
                invs.append(int(r[2]))
            #print len(null_ibd)
	    if self.pinv_sall==[]:
		self.pinv_sall=[0 for x in xrange(rep)]
                for idx,inv in enumerate(invs):
                    self.pinv_sall[inv]=null_ibd[idx][1]
	    if update_null:
		self.null_ibd = null_ibd
            #expectation and standard deviation
            ibd_pair=[x[0] for x in null_ibd]
            ibd_all=[x[1] for x in null_ibd]
            self.null_mean = float('%.9f'%(sum(ibd_pair)/len(ibd_pair)))
            self.null_std = self.std(ibd_pair,self.null_mean,False)
            self.sall_null_mean = float('%.9f'%(sum(ibd_all)/len(ibd_all)))
            self.sall_null_std = self.std(ibd_all,self.sall_null_mean,False)

        def null_perfect_rvibd(self,n_jobs,perfect_max,sall_flag=False,infer_flag=2):
        #null distribution of IBDs under perfect data
            screen_output.run_out("calculating theoretical null distribution for rvibd...")
            rep=pow(2,2*len(self.nonfounder))
            #print rep
            #Get all possible parental genotype configurations
            foundergt_list=[]
            foundergt_dic={}
            for f in self.founder:
                foundergt_list += self.fam_dict[f]['gt'][2*self.mid:2*self.mid+2]
            try:
                #if there are missing
                gt_id=self.conditional_prob['~combined'][1]
                for gt,prob in self.conditional_prob['~combined'][0].iteritems():
                    for f in set(self.missing_all)&set(self.founder):
                        f_idx=self.founder.index(f)
                        idx=gt_id.index(f)
                        foundergt_list[2*f_idx:2*f_idx+2]=[ord(x)-96 for x in gt[2*idx:2*idx+2]]
                    key_fgt=[]
                    for i in range(len(self.founder)):
                        key_fgt.extend(sorted(foundergt_list[2*i:2*i+2]))
                    key_fgt=tuple(key_fgt)
                    if key_fgt not in foundergt_dic:
                        foundergt_dic[key_fgt]=prob
                    else:
                        foundergt_dic[key_fgt]+=prob
            except:
                foundergt_dic[tuple(foundergt_list)]=1
            mean_pair, mean_all, var_pair, var_all = 0, 0, 0, 0
            all_ibd = manager.dict(self.all_ibd)
            if self.pinv_sall==[]:
                pinv_sall = Array(c_double,[0 for x in xrange(pow(2,2*len(self.invnf)))])
            else:
                pinv_sall=(c_double*len(self.pinv_sall))(*self.pinv_sall)
            pinv_key_dict = manager.dict(self.pinv_key_dict)
            pinv_pair_dict = manager.dict(self.pinv_pair_dict)
            combined_dist={}
            prob_and_mean=[]
            if self.simple:
                local_perfect_max=int(perfect_max/len(foundergt_dic.keys()))
            else:
                local_perfect_max=min(int(perfect_max/len(foundergt_dic.keys())),5)
            for fgt,prob in foundergt_dic.iteritems():
                #For each possible parental genotypes
                #Calculate corresponding IBD given inheritance vector
                #full_permutation=list(permutations(list(fgt),len(self.founder)*2))
                sample_flag=False
                full_permutation=[list(fgt)]
                if len(set(fgt))==1:
                    full_permutation_inv_raw=[(0,list(fgt))]
                else:
                    full_permutation_inv_raw=list(product(xrange(rep),full_permutation))
                ###parallel processing###
                inqueue = multiprocessing.Queue()
                null_ibd = manager.list([])
                if len(full_permutation_inv_raw) > local_perfect_max:
                    #sample_flag=True
                    full_permutation_inv = random.sample(full_permutation_inv_raw,local_perfect_max)
                else:
                    full_permutation_inv = full_permutation_inv_raw
                len_total=len(full_permutation_inv)
                for i in full_permutation_inv:
                    inqueue.put(i)
                #start n seperate processes
                procs = []
                for proc in range(n_jobs):
                    p = myProcess(proc,self.founder,self,inqueue,all_ibd,\
                        null_ibd,pinv_sall,pinv_key_dict,pinv_pair_dict,sall_flag,infer_flag)
                    p.start()
                    procs.append(p)
                    inqueue.put(None)
                TIMEOUT = 3600
                start = time.time()
                while time.time() - start <= TIMEOUT or len(null_ibd)<2:
                    if any(p.is_alive() for p in procs):
                        time.sleep(.1)  # Just to avoid hogging the CPU
                    else:
                        break
                else:
                    try:
                        #print("timed out, killing all processes")
                        for p in procs:
                            p.terminate()
                            p.join()
                    except:
                        pass
                while not inqueue.empty():
                    inqueue.get()
                tmp_ibd_pair = [x[0] for x in null_ibd]
                tmp_mean_pair = float('%.9f'%(sum(tmp_ibd_pair)/len_total))
                tmp_std_pair = self.std(tmp_ibd_pair,tmp_mean_pair,sample_flag)
                mean_pair += prob* tmp_mean_pair
                var_pair += prob*tmp_std_pair**2
                tmp_mean_all, tmp_std_all = 0, 0
                if sall_flag:
                    tmp_ibd_all = [x[1] for x in null_ibd]
                    tmp_mean_all = float('%.9f'%(sum(tmp_ibd_all)/len_total))
                    tmp_std_all = self.std(tmp_ibd_all,tmp_mean_all,sample_flag)
                    mean_all += prob* tmp_mean_all
                    var_all += prob*tmp_std_all**2
                prob_and_mean.append((prob,tmp_mean_pair,tmp_mean_all))
                dist_s=self.distribution(pall_flag=sall_flag,null_ibd=null_ibd)
                #if len(set(fgt))==1:
                #    print tmp_mean_pair, tmp_mean_all, tmp_std_pair, tmp_std_all, dist_s
                for v_idx,ibd_v in enumerate(dist_s[1]):
                    if ibd_v in combined_dist:
                        combined_dist[ibd_v]+=dist_s[0][v_idx]*prob
                    else:
                        combined_dist[ibd_v]=dist_s[0][v_idx]*prob
            ibd_keys = combined_dist.keys()
            self.dist_s = ([combined_dist[k] for k in ibd_keys], ibd_keys)
            self.all_ibd = all_ibd
            for tmp_ele in prob_and_mean:
                var_pair += tmp_ele[0]*(tmp_ele[1]-mean_pair)**2
                if sall_flag:
                    var_all += tmp_ele[0]*(tmp_ele[2]-mean_all)**2
            self.null_mean = mean_pair
            self.sall_null_mean = mean_all
            self.null_std = math.sqrt(var_pair)
            self.sall_null_std = math.sqrt(var_all)

	def stat(self):                         #get Z-score and observed IBD
		observe = self.ibd_total
		if self.null_mean:
		    try:
		        return (observe-self.null_mean)/self.null_std,observe #z score;observed IBD 
		    except ZeroDivisionError:
                        return 0, observe
		else:
		    return 0, observe
                
        def distribution(self,pall_flag,null_ibd=[]):             #get the IBD distribution under null hypothesis
            if null_ibd == []:
                null_ibd = self.null_ibd
            if null_ibd:
                if pall_flag:
                    null_ibd = [v for v in null_ibd if v[1]!=0]
                total = len(null_ibd)
                dis = []
                values = list(set(null_ibd))
                dis=[null_ibd.count(v)/total for v in values]
                return dis, values
            else:
                return [1],[self.ibd_total]        

def smart_float(x):
    try:
        return float('%.9f'%x)
    except:
        return None

def randomsample(dic,n_jobs,rep):
    #random sample in each family for overall IBD sharing
    rep=int(rep)
    famnum = len(dic)
    ibd_dist = []
    fam_struct = []      #unique family structure with null distribution
    structnum = []
    for i in range(famnum):
        if dic[i] not in fam_struct:
            fam_struct.append(dic[i])
            structnum.append(dic.count(dic[i]))
    ##multiple process####
    p = Pool(processes=n_jobs)
    #total_jobs=[(structnum,fam_struct) for x in xrange(rep)]
    result = p.imap(ransample,repeat((structnum,fam_struct),rep), chunksize=rep//n_jobs) 
    #result.wait()
    p.close()
    p.join()
    ibd_dist_raw = [(smart_float(r[0]),smart_float(r[1])) for r in result]
    ibd_dist = zip(*ibd_dist_raw)
    ibd_dist = [sorted(list(ibd_dist[0]),reverse=True),sorted(list(ibd_dist[1]),reverse=True)]
    return ibd_dist

def smart_int(x):
    try: 
        return int(x)
    except ValueError: 
        return 0

def likelihood_diff(x,*coeff):
    result=0
    c_list=[]
    try:
        c_list=list(coeff[0])
    except:
        c_list=list(coeff)
    for i in list(c_list):
        result+=i/(1+i*x)
    return result

def likelihood(coeff,x):
    result=0
    for i in list(coeff):
        result+=math.log(1+i*x)
    return result

def K_C(coeff,min_b):
	#Kong&Cox extension of asymptotic p-value
	x=Symbol('x',real=True)
	try:
            delta=list(solveset(likelihood_diff(x,tuple(coeff)),x,Interval(0,min_b)))
        except:
            #use numerical analysis to approach root if there's convergence problem by sympy
            delta_f=fsolve(likelihood_diff,0,args=tuple(coeff)) 
            delta=[tmp_delta for tmp_delta in list(delta_f) if tmp_delta>=0 and tmp_delta<=min_b]
	max_l=None
	#print "delta:"+repr(delta)
	if len(delta)==0:
	    #if there is no zero for likelihood_diff
	    #note likelihood_diff(0,tuple(coeff)) equals to original Z-score and should always be positive here.
	    try:
                d2=likelihood_diff(min_b,tuple(coeff))
		likelihood(tuple(coeff),min_b)
            except:
                drift=0.01*min_b
                min_b=min_b-drift
                d2=likelihood_diff(min_b,tuple(coeff))
            while d2<-0.1:
                min_b=0.99*min_b
                d2=likelihood_diff(min_b,tuple(coeff))
            max_l=likelihood(tuple(coeff),min_b)
        else:
            if len(delta)>1:
                l_candidates=[]
                for tmp_delta in delta:
                    try:
                        tmp_l=likelihood(tuple(coeff),tmp_delta)
                    except:
                        tmp_l=likelihood(tuple(coeff),0.99*tmp_delta)
                    l_candidates.append(tmp_l)
                max_l=max(l_candidates)
            else:
                try:
                    max_l=likelihood(tuple(coeff),delta[0])
                except:
                    max_l=likelihood(tuple(coeff),0.99*delta[0])
        #print max_l,likelihood(tuple(coeff),0)
        lrt=max_l-likelihood(tuple(coeff),0)
        if lrt<0:
            lrt=0
        lrt_stat=math.sqrt(2*lrt)
        return lrt_stat,stats.norm.sf(lrt_stat)

def execute(args):
	screen_output.log_out(args.path)
	famstruct = []                      #family structures
	famstruct_null = {}                 #famstruct_null={struct_id:[null_mean,null_std,null_ibd]}
	famstruct_perfect = []
	famstruct_null_perfect = {}
	sall_struct = {}
	for (dirpath,dirname,filename) in os.walk(args.path):
		fped = [f for f in filename if re.search(r'ped',f)]
		fdat = [d for d in filename if re.search(r'dat',d)]
		if len(fped) == 1 and len(fdat) == 1:
		    pedfile = fped[0]
		    markerfile = fdat[0]
		    tmp_dir=os.path.dirname(dirpath)
		    pdirname = os.path.dirname(tmp_dir)
		    prefix_name = os.path.basename(tmp_dir)
		    foldername="_".join(os.path.relpath(tmp_dir,args.path).split('/'))
		    if foldername == "":
			foldername=prefix_name
		    chr_m = re.search(r'(chr\w+)',markerfile)
		    chr_id = chr_m.group(1)                   #get the chromosome id
		    genename=os.path.basename(pdirname)#foldername.split('_gnomAD')[0]
		    screen_output.log_out("folder name:{}, chromosome:{}".format(foldername,chr_id))
		    #get the marker names
		    markername = []
		    markerfile = '{}/'.format(dirpath)+markerfile
		    pedfile = '{}/'.format(dirpath)+pedfile
		    cache_pdirname = '{}/cache/CACHE'.format(pdirname)
		    mfile = file(markerfile,'r').read()
		    mlines = mfile.split('\n')[1:-1]
		    for element in mlines:
			markername.append(element.split()[-1])
		    #output original result
		    if args.output is None:
			args.output = pdirname
		    result = file("{}/original_result.txt".format(args.output),'a')
		    rtext = '%s:\n'%genename
		    #take in ped files
		    f = file(pedfile,'r').read()
		    lines = f.split('\n')
		    lines.pop()
		    family = {}
		    fam_class = {}
		    fam_ids = []
		    marker_n=0
		    for line in lines:
			tmp=line.split()
			tmp[4:] = [smart_int(i) for i in tmp[4:]]        #get every element from line
			fid = tmp[0]                                   #family id
			if fid not in fam_ids:
			    fam_ids.append(fid)
			    family[fid]=[]
			marker_n = int((len(tmp)-6)/2)
			#family=[[[fam1member1][fam1member2]...],[fam2member1][]...]
			family[fid].append(tmp[1:])
		    fam_num = len(family)
		    #Extract marker freq info from cache
		    archive=ZipFile('{}/cache/{}.cache'.format(pdirname,prefix_name),'r')
		    if args.snv:
			maf_by_marker={}
			try:
			    with archive.open('CACHE/{}.{}.freq'.format(prefix_name,chr_id)) as cache_fh:
				for line in cache_fh:
				    tmp_elements=line.split()
				    mname=tmp_elements[1]
				    tmp_freq=[float(x) for x in tmp_elements[2:]]
				    if mname not in maf_by_marker:
					maf_by_marker[mname]=tmp_freq
			except:
			    pass
		    else:
			#get the freq file for CHP markers
			CHP_freq_byfam={}
			for mk_name in markername:
			    CHP_freq_byfam[mk_name]={}
			try:
			    # read from freq file to get freq information for CHP alleles
			    with archive.open('CACHE/{}.{}.freq'.format(prefix_name,chr_id)) as cache_fh:
				for line in cache_fh:
				    tmp=line.split()
				    CHP_freq=[]
				    tmp_freqs=[float(x) for x in tmp[2:]]
				    sum_freq=float('%.7f'%sum(tmp_freqs))
				    if len(tmp_freqs)>2 and tmp_freqs[0]==tmp_freqs[-1]:
					#when there is a missing allele denoted as '0' 
					#resulted from imputation step in SEQLinkage
					#use the freqs of alleles other than '0'
					tf_idxs=[tf_idx for tf_idx,t_freq in enumerate(tmp_freqs) if t_freq==tmp_freqs[0]]
					CHP_freq=[freq/sum(tmp_freqs[:tf_idxs[1]]) for freq in tmp_freqs[:tf_idxs[1]]]
				    else:
					if sum_freq == 1:
					    #when the sum of freqs is 1, stating that only one RV 
					    #locus present here
					    CHP_freq=tmp_freqs
					elif sum_freq < 1:
					    #when the sum is smaller than 1
					    #stating that there are more than one RV loci present in family
					    #adjust the freqs so that the sum is 1
					    CHP_freq=[freq/sum_freq for freq in tmp_freqs]
					else:
					    CHP_freq=[freq/sum(tmp_freqs[:-1]) for freq in tmp_freqs[:-1]]
				    CHP_freq_byfam[tmp[1]][tmp[0]]=CHP_freq
			except:
			    pass
		   ###NPL-pair calculation###################
		    z_dist = []              #distribution of NPL z-score per family
		    s_dist = []              #distribution of s-score per family
		    sall_dist = []
		    for m in range(marker_n):
			z_dist.append([])
			s_dist.append([])
			sall_dist.append([])
		    observe = [0]*marker_n
		    z_sum = [0]*marker_n
		    zall_sum = [0]*marker_n
		    z_pair = [[] for x in xrange(marker_n)]
		    z_all = [[] for x in xrange(marker_n)]
		    prior_b_pair = [[] for x in xrange(marker_n)]
		    prior_b_all = [[] for x in xrange(marker_n)]
		    asymp_pv = [1]*marker_n
		    Zpair_pr = [1]*marker_n
		    z_scores = [1]*marker_n
                    fam_conditional_prob=[]
                    fam_conditional_prob_null={}
                    fam_conditional_prob_perfect=[]
                    fam_conditional_prob_null_perfect={}
		    
		    markers_to_analyze=range(marker_n)
		    fam_to_analyze={}
		    for mkid in markers_to_analyze:
			fam_to_analyze[markername[mkid]]=copy.copy(fam_ids)
			fam_class[mkid]={}
		    for m in markers_to_analyze:
		        wt_fam_count=0
			screen_output.run_out("analyzing "+repr(markername[m]))
			if None in fam_to_analyze[markername[m]]:
			    continue
			for fid in fam_to_analyze[markername[m]]:
			    #for each family, construct a class of Family
			    print "fid:"+repr(fid)
			    fam = Family(m)
			    fam.rvibd = args.rvibd
			    if args.snv and maf_by_marker != {}:
				fam.mp_freq = maf_by_marker[markername[m]]
			    elif not args.snv:
				if CHP_freq_byfam[markername[m]].has_key(fid):
				    fam.mp_freq=CHP_freq_byfam[markername[m]][fid]
			    fam.setdict(family[fid])                   #set family dictionary
			    try:
				if args.info_only and not fam.info:
				    screen_output.run_out("uninformative family detected")
				    fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
				    fam.clean()
				    continue
				fam.classify_affect()
			    except:
				fam.err=True
			    if fam.err:
				screen_output.err_out("Error detected in fam {} for marker {}".format(fid,markername[m]))
				fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
				fam.clean()
				continue
			    fam_class[m][fid]=fam
			    if fam.wt_fam:
				wt_fam_count+=1
			    #print fam.info, fam.famstruct
			    #print fam.conditional_prob['~combined'][0]
			    if args.perfect:
                                if args.rvibd:
                                    try:
                                        fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                    except:
                                        fam_combined_info=(fam.famstruct,None)
                                    if fam.missing_all!=[]:
                                        if fam_combined_info not in fam_conditional_prob_perfect:
                                            fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=args.sall)
                                            fam_conditional_prob_perfect.append(fam_combined_info)
                                            fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
                                        else:
                                            fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std=fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]
                                    else:
                                        if fam.famstruct not in famstruct_perfect:
                                            fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=args.sall)
                                            famstruct_perfect.append(fam.famstruct)
                                            famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
                                        else:
                                            fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]
                                else:
                                    pfamstruct=copy.copy(fam.famstruct)
                                    pfamstruct.pop('info',None)
                                    if pfamstruct not in famstruct_perfect:
                                        fam.null_perfect(n_jobs=args.n_jobs,perfect_max=args.perfect_max)
                                        famstruct_perfect.append(pfamstruct)
                                        famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]=[fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]
                                    else:
                                        fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]
			    else:
				#use permutations to calculate nullibd
				if fam.missing_all!=[] and not fam.wt_fam:
                                    fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                    if fam_combined_info not in fam_conditional_prob:
                                        fam.nullibd(rep=50,n_jobs=args.n_jobs,sall_flag=args.sall)
                                        fam_conditional_prob.append(fam_combined_info)
                                        fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=[fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]
                                    else:
                                        fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std=fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]
				elif fam.info and not fam.wt_fam:   
				    if fam.famstruct not in famstruct:
                                        fam.nullibd(rep=100,n_jobs=args.n_jobs,sall_flag=args.sall)
                                        famstruct.append(fam.famstruct)
                                        fam.dist_s=fam.distribution(pall_flag=args.sall)
                                        famstruct_null[str(famstruct.index(fam.famstruct))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]
                                    else:
                                        fam.dist_s,fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std=famstruct_null[str(famstruct.index(fam.famstruct))]
			    ###############
			    allele_key=''
			    for person in sorted(fam.fam_dict.keys()):
				if person in fam.missing_all:
				    allele_key+='00'
				else:
				    allele_key+=''.join(map(str,fam.fam_dict[person]['gt'][2*fam.mid:2*fam.mid+2]))
			    if allele_key in fam.all_ibd:
				fam.ibd_total, fam.ibd_sall=fam.all_ibd[allele_key]
			    else:
				try:
				    fam.ibd_total,fam.ibd_sall=fam.cal_ibd(sall_flag=args.sall)
				except:
				    fam.err=True
			    if fam.err:
				screen_output.err_out("Error detected in fam {} for marker {}".format(fid,markername[m]))
				fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
				fam.clean()
				continue
			    z,o = fam.stat()                             #get Z-score and observed IBD among affected individuals for each family 
			    if not fam.wt_fam:
				z_sum[m] += z
			    z_pair[m].append(z)
			    #print fam.ibd_total
			    #print "z:"+repr(z)
			    observe[m] = float('%.9f'%(observe[m]+o))                     #combined value of IBD across family
                            if args.sall:
                            ###S-all statistic########
				z_sall=0
				try:
				    z_sall=(fam.ibd_sall-fam.sall_null_mean)/fam.sall_null_std
				except ZeroDivisionError:
				    z_sall=0
				#print fam.ibd_sall,fam.sall_null_mean,fam.sall_null_std
				if not fam.wt_fam:
				    zall_sum[m] += z_sall
				#print "zall:"+repr(z_sall)
				z_all[m].append(z_sall)
                                ###############
                            if args.perfect:
                                if args.kc:
                                    b_pair=fam.null_std/(fam.null_mean-min([x[0] for x in fam.null_ibd]))
                                    prior_b_pair[m].append(b_pair)
                                    if args.sall:
                                        min_sall=min([x[1] for x in fam.null_ibd])
                                        b_all=fam.sall_null_std/(fam.sall_null_mean-min_sall)
                                        prior_b_all[m].append(b_all)
			#asymptotic p-value NPL-pair
			#NPL-pair score; weighted by sqrt of family number
			fam_num = len([x for x in fam_to_analyze[markername[m]] if x != None])-wt_fam_count
			if fam_num == 0:
			    asymp_pv[m] = (0.5,0.5)
			    z_scores[m] = (0,0)
			    continue
			Z_pair = z_sum[m]/math.sqrt(fam_num)
			p_norm = stats.norm.sf(Z_pair)    #asymptotic p-value
			#print Z_pair,p_norm
			asymp_pv[m] = (p_norm,p_norm)
			z_scores[m] = (Z_pair,Z_pair)
			if args.sall:
			    #asymptotic p-value NPL-all
			    Z_all = zall_sum[m]/math.sqrt(fam_num)
			    pall_norm = stats.norm.sf(Z_all)
			    #print Z_all,pall_norm
			    asymp_pv[m] = (p_norm,pall_norm)
			    z_scores[m] = (Z_pair,Z_all)
			if args.kc and z_scores[m][0]>0 and z_scores[m][1]>0:
			    ####Kong&Cox#####
			    ##only tested for delta>0###
			    #print fam_num
			    coeff_pair=[]
			    for tmp_z_pair in z_pair[m]:
				coeff_pair.append(tmp_z_pair/math.sqrt(fam_num))
			    min_b_pair=min(prior_b_pair[m])*math.sqrt(fam_num)
			    #print "min_b_pair:"+repr(min_b_pair)
			    Zpair_kc,p_pair_norm_kc=K_C(coeff_pair,min_b_pair)	
			    #print Zpair_kc,p_pair_norm_kc
			    asymp_pv[m] = (p_pair_norm_kc,p_pair_norm_kc)
			    z_scores[m] = (Zpair_kc,Zpair_kc)
			    if args.sall:
				coeff_all=[]
				for tmp_z_all in z_all[m]:
				    coeff_all.append(tmp_z_all/math.sqrt(fam_num))
				min_b_all=min(prior_b_all[m])*math.sqrt(fam_num)
				#print "min_b_all:"+repr(min_b_all)	
				Zall_kc,p_all_norm_kc=K_C(coeff_all,min_b_all)	
				#print Zall_kc,p_all_norm_kc
				asymp_pv[m] = (p_pair_norm_kc,p_all_norm_kc)
				z_scores[m] = (Zpair_kc,Zall_kc)
		    pr_s_min, pr_sall_min = 1, 1
		    small_num=3
		    if min([x for x in asymp_pv if x !=1])[0]>0.5:
			small_num=1
		    for tapv in sorted([(x,idx_apv) for idx_apv,x in enumerate(asymp_pv) if x !=1])[:small_num]:
			#take the 3 smallest asymptotic p-value to calculate empiric p-value
			apv=tapv[0]
			idx_apv = tapv[1]
			asymp_z_pair,asymp_z_all=z_scores[idx_apv]
			rtext += 'marker: %s\n'%markername[idx_apv]
			rtext += 'asymptotic p-value:'+repr(apv)+'\n'
			rtext += 'Z-pair:{}\tZ-all:{}\n'.format(asymp_z_pair,asymp_z_all)
			fam_num = len([x for x in fam_to_analyze[markername[idx_apv]] if x != None])
			pr_s, pr_sall = None, None
			print "marker:%s"%markername[idx_apv]
			print apv
			rep=args.rep
			fam_rep=args.fam_rep
			fam_reps=[50, 100, 1000, args.fam_rep]
                	reps=[500, 5000, 50000, args.rep]
			if fam_num==0:
			    pr_s, pr_sall = 0.5, 0.5
			if apv[0]>=args.cut and apv[1]>=args.cut:
			    #Adaptive permutations
			    #Use asymptotic p-values if the results are not pointwise significant
			    pr_s, pr_sall = apv[0],apv[1]
			elif args.exact and not args.force:
			    if apv[0]>0.5:
				fam_rep=fam_reps[0]
				rep=reps[0]
			    elif apv[0]>0.1:
				fam_rep=fam_reps[1]
				rep=reps[1]
			    elif apv[0]>0.01:
				fam_rep=fam_reps[2]
				rep=reps[2]
			while (pr_s,pr_sall)==(None,None) or pr_s==1/(1+rep) or pr_sall==1/(1+rep):
			    #permutations needed for empirical p-values
			    if (pr_s,pr_sall)==(None,None):
				pr_s, pr_sall = 1, 1
			    cz_dist = []        #combined distribution of Z-score of S-pair statistic
			    csall_dist = []     #combined distribution of Z-score of S-all statistic
			    z_dist[idx_apv]=[]
			    sall_dist[idx_apv]=[]
			    p_norm,pall_norm = apv[0],apv[1]
			    z_sum_precise=0
			    zall_sum_precise=0
			    increase_permutation=False
			    if pr_s==1/(1+rep) or pr_sall==1/(1+rep):
				#if empiric p-value reached limit, more permutation is needed
				if fam_rep==args.fam_rep:
				    break	
				increase_permutation=True
				#print "increase permutation"
				#print fam_rep,rep
				idx=fam_reps.index(fam_rep)
				fam_rep=fam_reps[idx+1]
				rep=reps[idx+1]
				#print fam_rep,rep
			    pall_flag=True if fam_rep>1000 or args.sall else False
			    infer_full_flag=1 if fam_rep>1000 else 2
			    rtext+='permutations:{}\n'.format(rep)
			    for fid in fam_to_analyze[markername[idx_apv]]:
				if fid != None:
				    fam = fam_class[idx_apv][fid]
				    if not increase_permutation:
					fam.null_ibd = []
				    if not fam.simple:
					infer_full_flag=2
				    if infer_full_flag==1 and fam.missing_all!=[] and not fam.wt_fam:
					fam.GT_infer(infer_flag=infer_full_flag)
					fam.ibd_total,fam.ibd_sall=fam.cal_ibd(sall_flag=args.sall)
					fam.all_ibd={}
				    if fam.missing_all!=[]:
                                        fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                        if args.perfect and not args.rvibd:
                                            if fam_combined_info not in fam_conditional_prob:
                                                fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True)
                                                fam_conditional_prob.append(fam_combined_info)
                                                fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=fam.null_ibd
                                            else:
                                                fam.null_ibd=fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]
                                                if pall_flag and fam.null_ibd[-1][1]==0 or len(fam.null_ibd)<fam_rep:
                                                    fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True)
                                                    fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=fam.null_ibd
					elif args.rvibd:
					    if args.perfect:
						fam.null_ibd=[]
						try:
						    fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std=fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]
						except:
                                            	    fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag)
						    fam_conditional_prob_perfect.append(fam_combined_info)
						    fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
						if pall_flag and fam.sall_null_mean==0:
                                            	    fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag)
						    fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
					    elif fam.null_std==0:
						#unlikely to get different values through permutations
						#this can only occur in rvibd when missing
						fam.null_ibd=[]
						if fam_combined_info not in fam_conditional_prob_perfect:
                                            	    fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag)
						    fam_conditional_prob_perfect.append(fam_combined_info)
						    fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
						else:
						    fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std=fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]
						    if pall_flag and fam.sall_null_mean==0:
                                            	        fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag)
							fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
                                            else: #use permutations
						fam.dist_s=None
                                                if fam_combined_info not in fam_conditional_prob:
                                                    fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag)
                                                    fam_conditional_prob.append(fam_combined_info)
                                                    fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=[fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]
                                                else:
                                                    fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std=fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]
                                                    if pall_flag and fam.null_ibd[-1][1]==0 or len(fam.null_ibd)<fam_rep:
                                                        fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag)
                                                        fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=[fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]
                                        if fam.null_ibd:
                                            add_time=0
                                            while max([x[0] for x in fam.null_ibd])<fam.ibd_total:
                                                #print "add additional permutations"
                                                if args.perfect and not args.rvibd:
                                                    fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True)
                                                else:
                                                    fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag)
                                                add_time+=1
                                                if add_time > 10:
                                                    break
                                            #print fam.ibd_total,sorted(fam.null_ibd)[-10:]
                                            #print len(fam.null_ibd)
				    elif fam.info:
					if args.perfect and not args.rvibd:
					    if fam.famstruct not in famstruct:
						fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True)
						famstruct.append(fam.famstruct)
						fam.dist_s=fam.distribution(pall_flag=pall_flag)
                                                famstruct_null[str(famstruct.index(fam.famstruct))]=[fam.dist_s,len(fam.null_ibd)]
					    else:
                                                fam.dist_s,fam.null_ibd_len=famstruct_null[str(famstruct.index(fam.famstruct))]
                                                if pall_flag and fam.dist_s[1][0][1]==0 or fam.null_ibd_len<fam_rep:
                                                    fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True)
                                                    fam.dist_s=fam.distribution(pall_flag=pall_flag)
                                                    famstruct_null[str(famstruct.index(fam.famstruct))]=[fam.dist_s,len(fam.null_ibd)]
					else:
                                            if args.perfect:
                                                fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]
						if pall_flag and fam.sall_null_mean==0:
                                            	    fam.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag)
						    famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std] 
					    else:
                                                if fam.famstruct not in famstruct:
                                                    fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag)
                                                    famstruct.append(fam.famstruct)
                                                    fam.dist_s=fam.distribution(pall_flag=pall_flag)
                                                    famstruct_null[str(famstruct.index(fam.famstruct))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]
                                                else:
                                                    fam.dist_s,fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std=famstruct_null[str(famstruct.index(fam.famstruct))]
                                                    if pall_flag and fam.dist_s[1][0][1]==0 or len(fam.null_ibd)<fam_rep:
                                                        fam.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag)
                                                        fam.dist_s=fam.distribution(pall_flag=pall_flag)
                                                        famstruct_null[str(famstruct.index(fam.famstruct))]=[fam.dist_s,fam.null_mean,fam.null_std,fam.null_ibd,fam.sall_null_mean,fam.sall_null_std]				
				    if args.rvibd:
					pfamstruct=copy.copy(fam.famstruct)
					pfamstruct.pop('info',None)
					if pfamstruct not in famstruct_perfect:
                                            fam.null_perfect(n_jobs=args.n_jobs,perfect_max=args.perfect_max,update_null=False)
					    famstruct_perfect.append(pfamstruct)
					    famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]=[fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std]
					else:
					    fam.null_mean,fam.null_std,fam.sall_null_mean,fam.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]
				    if pall_flag:
					if fam.ibd_sall == 0:
					###S-all statistic########
					    fam.ibd_sall=fam.cal_ibd({},fam.conditional_prob,sall_flag=pall_flag)[1]	
				    tmp = {}
				    if args.exact and cz_dist == []:
					d,v = fam.dist_s if fam.dist_s else fam.distribution(pall_flag=pall_flag)            #get the expected IBD distribution
					#print fid
					#print [d[v.index(tv)] for tv in sorted(v)],sorted(v)
					#print sum([d[v.index(tv)] for tv in sorted(v,reverse=True) if tv[0]>=fam.ibd_total])
					#print fam.ibd_total, fam.null_mean, fam.null_std
					if fam.null_mean:
					    try:
						zv_pair = (fam.ibd_total-fam.null_mean)/fam.null_std
					    except ZeroDivisionError:
						zv_pair = 0
					else:
					    zv_pair = 0
					z_sum_precise+=float('%.9f'%zv_pair)
					if pall_flag:
					    if fam.sall_null_mean:
						try:
						    zv_all = (fam.ibd_sall-fam.sall_null_mean)/fam.sall_null_std
						except ZeroDivisionError:
						    zv_all = 0
					    else:
						zv_all = 0
					    zall_sum_precise+=float('%.9f'%zv_all)
					for i in range(len(v)):
					    ibd_pair, ibd_all = v[i]
					    if fam.null_mean:
						try:
						    zv_pair = (ibd_pair-fam.null_mean)/fam.null_std
						except ZeroDivisionError:
						    zv_pair = 0
					    else:
						zv_pair = 0
					    zv_pair=float('%.9f'%zv_pair)
					    if pall_flag:
						if fam.sall_null_mean:
						    try:
							zv_all = (ibd_all-fam.sall_null_mean)/fam.sall_null_std
						    except ZeroDivisionError:
							zv_all = 0
						else:
						    zv_all = 0
						zv_all=float('%.9f'%zv_all)
					    else:
						zv_all=None
					    if (zv_pair,zv_all) not in tmp:
						tmp[(zv_pair,zv_all)]=0
					    tmp[(zv_pair,zv_all)] += d[i]                  #{z score: probability}
					z_dist[idx_apv].append(tmp)
		    ###################################################
			    if cz_dist == []:
				if args.exact:              #random sample
                                    screen_output.run_out('randomly sampling..')
                                    cz_dist,csall_dist=randomsample(z_dist[idx_apv],args.n_jobs,rep)
			    if args.exact:
				z_sum_precise=float('%.9f'%z_sum_precise)
				for i_dx,i in enumerate(cz_dist):
				    if i >= z_sum_precise:
					continue
				    else:
					pr_s=(i_dx+1)/(len(cz_dist)+1)
					break
			    #NPL-pair score; weighted by sqrt of family number
			    Z_pair = z_sum_precise/math.sqrt(fam_num)                  
			    lod_pair = Z_pair**2/(2*math.log(10))
			    rtext += 'S-pair:\n'
			    rtext += 'asymptotic p-value:'+repr(apv[0])+'\n'
			    if args.exact:
				rtext += 'empiric p-value:'+repr(pr_s)+'\n'
			    #############S-all#######################
                            if args.exact and None not in csall_dist:
				zall_sum_precise=float('%.9f'%zall_sum_precise)
				for i_dx,i in enumerate(csall_dist):
				    if i >= zall_sum_precise:
					continue
				    else:
					pr_sall=(i_dx+1)/(len(csall_dist)+1)
					break
			    else:
				pr_sall = pr_s
			    rtext += "S-all:\n"
			    if zall_sum_precise:
				Z_all = zall_sum_precise/math.sqrt(fam_num)
				#pall_norm = stats.norm.sf(Z_all)
			    if args.exact:
				rtext += "empiric p-value:"+repr(pr_sall)+"\n"
			    rtext += "asymptotic p-value:"+repr(apv[1])+"\n"
			    ####################
			pr_s_min = pr_s if pr_s<pr_s_min else pr_s_min
			pr_sall_min = pr_sall if pr_sall<pr_sall_min else pr_sall_min
			possible_min = 1/(int(rep)+1)
			if (pr_s_min,pr_sall_min) == (possible_min,possible_min):
			    break
		    result.write(rtext)
		    result.close()
		    pv = file("{}/pvalue.txt".format(args.output),'a')
		    marker_num=len(markername)
		    if args.exact:
			pv.write('%s\t%.9f\t%.9f\n'%(genename,pr_s_min,pr_sall_min))
		    else:
			p_pair_min,p_all_min=min([x for x in asymp_pv if x!=1])
			pv.write('%s\t%.9f\t%.9f\n'%(genename,p_pair_min,p_all_min))
		    pv.close()

