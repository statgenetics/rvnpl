#!/usr/bin/python 
from __future__ import division
import os
import re
import random
import copy
from zipfile import ZipFile
from ctypes import *
from itertools import repeat,product,permutations
from multiprocessing import Process, Value, Manager, Pool
import multiprocessing
import math
import numpy as np
import time
import string
from scipy import stats
from RVNPLcpp import cpostInv
from RVNPL import screen_output
from RVNPL.npl_class import NPL

class QNPL(NPL):
    def __init__(self,family):
        NPL.__init__(self,family)
        self.cond_null_ibd = {}
        #matrix 
        self.S = np.matrix([])
        self.D = np.matrix([])
        self.ES = np.matrix([])
        self.ED = np.matrix([])
        self.ET = np.matrix([])
        self.covSS = np.matrix([])
        self.covDD = np.matrix([])
        self.covSD = np.matrix([])
        self.covTT = np.matrix([])
        self.covY_inverse = np.matrix([])
        self.B = np.matrix([])
        self.H = np.matrix([])
        self.numerator = 0
        self.denominator = 0
    def cal_ibd_inv(self,v,fgt=[]):
        #calculate IBD among nonfounders under a specific inheritance vector
        #assuming fully informative markers
        t_alleles=self.family.assign_allele(v)
        shared_alleles=[]
        for t_pair in self.family.pairs:
            t_shared_alleles=0
            for tt_a in t_alleles[t_pair[0]]:
                if tt_a in t_alleles[t_pair[1]]:
                    if self.family.rvibd and tt_a in fgt or not self.family.rvibd:
                        t_shared_alleles+=1
            shared_alleles.append(t_shared_alleles/2)
        return shared_alleles         #pairwise IBD

    def combine_inherit_vec(self,m=-1,alleles={},conditional_prob={}):
        #calculate the overall posterior inv probability
        combined_postinv_prob = {}
        pairnum=len(self.family.pairs)
        self.family.post = [[0 for y in range(pairnum)] for x in range(pairnum)]
        ibdall=[0]*len(self.family.pairs)
        if m != -1:
            for iid in self.family.fam_dict.keys():
                alleles[iid] = self.family.fam_dict[iid]['gt'][2*m:2*m+2]
            conditional_prob=self.family.conditional_prob
        try:
            tmp_dic=conditional_prob['~combined'][0]
            gt_id=conditional_prob['~combined'][1]
        except:
            tmp_dic={'00':1}
        for tmp_gt,prob in tmp_dic.iteritems():
            # reassign genotypes
            for m_person in self.family.sorted_miss_persons:
                idx=gt_id.index(m_person[0])
                alleles[m_person[0]]=[ord(x)-96 for x in tmp_gt[2*idx:2*idx+2]]
            fgt,nfgt=[],[]
            if self.family.rvibd:
                for iid in self.family.fam_dict.keys():
                    if iid in self.family.founder:
                        fgt.extend(['{}{}'.format(iid,a_idx) for a_idx,tmp_a in enumerate(alleles[iid]) if tmp_a!=1])
                    else:
                        nfgt+=alleles[iid]
            if self.family.rvibd and (fgt==[] or set(nfgt)==set([1])):
                continue
            tmp_postinv=cpostInv.apply(self.family.invnf, self.family.founder, self.family.parents, alleles, self.family.rvibd)
            len_postinv = len(tmp_postinv)
            ibd_pair_all={}
            for tpostinv in set(tmp_postinv):
                ibd_pair=tuple(self.cal_ibd_inv(tpostinv,fgt))
                #print tpostinv, ibd_pair
                count=tmp_postinv.count(tpostinv)
                if ibd_pair not in ibd_pair_all:
                    ibd_pair_all[tuple(ibd_pair)]=count
                else:
                    ibd_pair_all[tuple(ibd_pair)]+=count
            #print "total:",len_postinv
            for tmp_pair_ibd in ibd_pair_all:
                count=ibd_pair_all[tmp_pair_ibd]
                tmp_prob=1/len_postinv*count*prob
                for idx1 in range(pairnum):
                    if tmp_pair_ibd[idx1]!=0:
                        ibdall[idx1]+=tmp_pair_ibd[idx1]*tmp_prob
                        for idx2 in range(idx1,pairnum):
                            self.family.post[idx1][idx2]+=tmp_pair_ibd[idx1]*tmp_pair_ibd[idx2]*tmp_prob
                            self.family.post[idx2][idx1]=self.family.post[idx1][idx2]
        if m != -1:
            self.family.ibdall = ibdall
        else:
            return ibdall

    def null_inv(self):
        ##get the pairwise IBD under each possible inheritance vector
        nfnum = len(self.family.nonfounder)
        pairnum = len(self.family.pairs)
        for i in xrange(pow(2,2*nfnum)):
            tmp = bin(i)[2:]
            tmpinv = '0'*(2*nfnum-len(tmp))+tmp
            self.family.prioribd[i]=self.cal_ibd_inv(tmpinv)      #pairwise IBD under each inv
        #print len(self.family.prioribd) 
        #########################
            for idx1 in range(pairnum):
                for idx2 in range(idx1,pairnum):
                    self.family.prior[idx1][idx2] += self.family.prioribd[i][idx1]*self.family.prioribd[i][idx2]/(2**(2*nfnum))
                    self.family.prior[idx2][idx1] = self.family.prior[idx1][idx2]

    def null_permute(self,n_jobs,perfect_max,conditional_prob=[]):
    ##get pairwise IBD under each prior inheritance vector while maintain founder genotype
        screen_output.run_out("calculating theoretical null distribution for rvibd...")
        paircount=len(self.family.pairs)
        self.family.expect_pair_ibd=[0 for x in range(paircount)]
        self.family.prior = [[0 for y in range(paircount)] for x in range(paircount)]
        rep=pow(2,2*len(self.family.nonfounder))
        #print rep
        #Get all possible parental genotype configurations
        foundergt_list=[]
        foundergt_dic={}
        if conditional_prob==[]:
            conditional_prob=self.family.conditional_prob
        for f in self.family.founder:
            foundergt_list += self.family.fam_dict[f]['gt'][2*self.family.mid:2*self.family.mid+2]
        try:
            #if there are missing
            gt_id=conditional_prob['~combined'][1]
            for gt,prob in conditional_prob['~combined'][0].iteritems():
                for f in set(self.family.missing_all)&set(self.family.founder):
                    f_idx=self.family.founder.index(f)
                    idx=gt_id.index(f)
                    foundergt_list[2*f_idx:2*f_idx+2]=[ord(x)-96 for x in gt[2*idx:2*idx+2]]
                key_fgt=[]
                for i in range(len(self.family.founder)):
                    key_fgt.extend(sorted(foundergt_list[2*i:2*i+2]))
                key_fgt=tuple(key_fgt)
                if key_fgt not in foundergt_dic:
                    foundergt_dic[key_fgt]=prob
                else:
                    foundergt_dic[key_fgt]+=prob
        except:
            foundergt_dic[tuple(foundergt_list)]=1
        local_perfect_max=int(perfect_max/len(foundergt_dic.keys()))
        p = Pool(processes=n_jobs)
        for fgt,prob in foundergt_dic.iteritems():
            #For each possible parental genotypes
            #Calculate corresponding IBD given inheritance vector
            sample_flag=False
            full_permutation=[list(fgt)]
            if len(set(fgt))==1:
                full_permutation_inv_raw=[(0,list(fgt))]
            else:
                full_permutation_inv_raw=xrange(rep) #list(product(xrange(rep),full_permutation))
            if len(full_permutation_inv_raw) > local_perfect_max:
                #sample_flag=True
                full_permutation_inv_tmp = random.sample(xrange(rep),local_perfect_max)
            else:
                full_permutation_inv_tmp = full_permutation_inv_raw
            if isinstance(full_permutation_inv_tmp[0],tuple):
                full_permutation_inv=full_permutation_inv_tmp
            else:
                full_permutation_inv=list(product(full_permutation_inv_tmp,full_permutation))
            len_total=len(full_permutation_inv)
            csize=len_total//n_jobs
            if csize<1:
                csize=1
            #p = Pool(processes=n_jobs)
            result = p.imap(null_generator_pairibd,[(i,self) for i in full_permutation_inv], chunksize=csize)
            pair_ibd=[r for r in result]
            #start n seperate processes
            for tmp_pair_ibd in set(pair_ibd):
                count=pair_ibd.count(tmp_pair_ibd)
                #print tmp_pair_ibd
                #print count
                tmp_prob=1/len_total*count*prob
                for pair_idx in range(paircount):
                    t_ibd=tmp_pair_ibd[pair_idx]*tmp_prob
                    if t_ibd !=0:
                        self.family.expect_pair_ibd[pair_idx]+=t_ibd
                        for pair_jdx in range(pair_idx,paircount):
                            self.family.prior[pair_idx][pair_jdx]+=t_ibd*tmp_pair_ibd[pair_jdx]
                            self.family.prior[pair_jdx][pair_idx]=self.family.prior[pair_idx][pair_jdx]
        p.close()
        p.join()

    def setmatrix(self):
        ##set the values for the matrix S, D and IBDm
        s = []
        es = []  #expectation of S
        d = []
        ed = []  #expectation of D
        tau = self.family.expect_pair_ibd if self.family.rvibd else []
        for aid,bid in self.family.pairs:
            #for each pair of relatives (include parent-offspring)
            trait1 = self.family.fam_dict[aid]['trait']
            trait2 = self.family.fam_dict[bid]['trait']
            s.append((trait1+trait2)**2)        #squared sum of traits
            d.append((trait1-trait2)**2)        #squared difference of traits
            if not self.ES.size:
                r = self.family.corre([aid,bid])
                es.append(2*(1+r))
                ed.append(2*(1-r))
                if not self.family.rvibd:
                    tau.append(self.family.expect([aid,bid]))  #expected IBD between a relative pair
        self.S = np.matrix(s)
        self.D = np.matrix(d)
        if not self.ES.size:
            self.ES = np.matrix(es)
            self.ED = np.matrix(ed)
            self.ET = np.matrix(tau)

    def setcov(self):
        ##set the covariance matrix
        pairnum=len(self.family.pairs)
        ss = np.zeros((pairnum,pairnum))
        dd = np.zeros((pairnum,pairnum))
        sd = np.zeros((pairnum,pairnum))
        tt = np.zeros((pairnum,pairnum))
        for idx1 in range(pairnum):
            for idx2 in range(idx1,pairnum):
                ##get the covariance matrix for tau(i,j),tau(k,l)
                prior,post = 0,0
                prior_expect = 0
                pair1 = self.family.pairs[idx1]
                pair2 = self.family.pairs[idx2]
                #print pair1, pair2
                if not self.family.rvibd:
                    prior_expect=self.family.expect(pair1)*self.family.expect(pair2)
                else:
                    prior_expect=self.family.expect_pair_ibd[idx1]*self.family.expect_pair_ibd[idx2]
                prior = self.family.prior[idx1][idx2]
                post = self.family.post[idx1][idx2]
                #print prior, post, self.ibdall[idx1], self.ibdall[idx2]
                #print self.expect(pair1), self.expect(pair2)
                #covariance for IBD of two pairs
                tt[idx2,idx1]=tt[idx1,idx2] = (prior-prior_expect)-\
                        (post-self.family.ibdall[idx1]*self.family.ibdall[idx2])
                if not self.covSS.size:
                    rik = self.family.corre([pair1[0],pair2[0]])
                    ril = self.family.corre([pair1[0],pair2[1]])
                    rjk = self.family.corre([pair1[1],pair2[0]])
                    rjl = self.family.corre([pair1[1],pair2[1]])
                    ss[idx2,idx1]=ss[idx1,idx2] = 2*(rik+ril+rjk+rjl)**2
                    #print pair1,pair2,rik,ril,rjk,rjl
                    dd[idx2,idx1]=dd[idx1,idx2] = 2*(rik+rjl-ril-rjk)**2
                    sd[idx1,idx2] = 2*(rik+rjk-ril-rjl)**2
                    sd[idx2,idx1] = 2*(rik+ril-rjk-rjl)**2
        self.covTT = np.matrix([tmt for tmt in tt])
        #print self.ibdall
        #print self.covTT
        #print self.ibdall
        if not self.covSS.size:
            self.covSS = np.matrix(ss)
            self.covDD = np.matrix(dd)
            self.covSD = np.matrix(sd)

    def reg(self):
        ##Regression Step
        pairnum = len(self.family.pairs)
        all_members=len(self.family.all_members)#int((math.sqrt(8*pairnum+1)+1)/2)
        #print pairnum, all_members
        Y= np.hstack([self.S,self.D[0,0:all_members]]).T
        Yc = Y-np.hstack([self.ES,self.ED[0,0:all_members]]).T
        Tc = (self.family.ibdall-self.ET).T
        if not self.covY_inverse.size:
            covss = self.covSS
            covsd = self.covSD[:,0:all_members]
            covdd = self.covDD[0:all_members,0:all_members]
            ###########
            covY = np.vstack([np.hstack([covss,covsd]), np.hstack([covsd.T,covdd])])
            self.covY_inverse=covY.getI()
        if not self.B.size:
            diag1 = np.diag([2]*pairnum)
            diag2 = np.diag([-2]*pairnum)
            self.H = np.hstack([diag1,diag2[:,0:all_members]])
            self.B = self.H*self.covY_inverse*Yc
        numerator = self.B.T*Tc
        denominator = self.B.T*self.covTT*self.B
        self.numerator = numerator[0,0]
        self.denominator = denominator[0,0]

    def matrix_cal(self,alleles,conditional_prob):
        #the whole matrix calculation step
        ##set the values for the matrix S, D and IBDm
        ###########################
        ibdall=self.combine_inherit_vec(alleles=alleles,conditional_prob=conditional_prob)
        #print ibdall
        ##set the covariance matrix
        pairnum=len(self.family.pairs)
        all_member_num=len(self.family.all_members)
        tt = np.zeros((pairnum,pairnum))
        for idx1 in range(pairnum):
            for idx2 in range(idx1,pairnum):
                ##get the covariance matrix for tau(i,j),tau(k,l)
                prior,post = 0,0
                prior_expect = 0
                pair1 = self.family.pairs[idx1]
                pair2 = self.family.pairs[idx2]
                #print pair1,pair2
                if not self.family.rvibd:
                    prior_expect=self.family.expect(pair1)*self.family.expect(pair2)
                else:
                    prior_expect=self.family.expect_pair_ibd[idx1]*self.family.expect_pair_ibd[idx2]
                prior = self.family.prior[idx1][idx2]
                post = self.family.post[idx1][idx2]
                tt[idx2,idx1]=tt[idx1,idx2] = (prior-prior_expect)-(post-ibdall[idx1]*ibdall[idx2])
                #print tt[idx1,idx2]
                #print [idx1,idx2]
        covTT = np.matrix(tt)
        #print covTT
        if self.family.rvibd:
            self.ET = np.matrix(self.family.expect_pair_ibd)
        Tc = (ibdall-self.ET).T
        B = self.B
        numerator = B.T*Tc
        denominator = B.T*covTT*B
        return numerator[0,0], denominator[0,0]

    def matrix_cal_trait(self,new_traits):
        #the whole matrix calculation step
        ##set the values for the matrix S, D and IBDm
        ###########################
        ##set the values for the matrix S, D and IBDm
        all_members=len(self.family.all_members)
        s = []
        d = []
        for aid,bid in self.family.pairs:
            #for each pair of relatives (include parent-offspring)
            trait1 = new_traits[aid]
            trait2 = new_traits[bid]
            s.append((trait1+trait2)**2)        #squared sum of traits
            d.append((trait1-trait2)**2)        #squared difference of traits
        S = np.matrix(s)
        D = np.matrix(d)
        Y= np.hstack([S,D[0,0:all_members]]).T
        Yc = Y-np.hstack([self.ES,self.ED[0,0:all_members]]).T
        B = self.H*self.covY_inverse*Yc
        Tc = (self.family.ibdall-self.ET).T
        numerator = B.T*Tc
        denominator = B.T*self.covTT*B
        return numerator[0,0], denominator[0,0]

    def nulldist(self,n_jobs,perfect_max,rep=5000):
        print "calculating null dist for a single fam.."
        nfnum = len(self.family.nonfounder)
        null_nd={}
        ###parallel processing###
        inqueue = multiprocessing.Queue()
        manager=Manager()
        null_ibd = manager.list([])
        all_ibd = manager.dict({})
        cond_null_ibd = manager.dict(self.cond_null_ibd)
        if self.conditional_prob:
            cond_key=self.conditional_prob['~combined'][0].items()[0]
            if cond_key not in cond_null_ibd.keys():
                cond_null_ibd[cond_key]=[self.expect_pair_ibd,self.prior]
        #print cond_null_ibd.keys()
        for i in xrange(rep):
            inqueue.put(i)
        #start n seperate processes
        procs = []
        for proc in range(n_jobs):
            p = myProcess(proc,self.family.founder,self,inqueue,all_ibd,cond_null_ibd,null_ibd,n_jobs,perfect_max)
            p.start()
            procs.append(p)
            inqueue.put(None)
        TIMEOUT = 3600
        start = time.time()
        while time.time() - start <= TIMEOUT:
            #print time.time()-start, len(null_ibd)
            if any([p.is_alive() for p in procs]) and len(null_ibd)<rep:
                time.sleep(.1)  # Just to avoid hogging the CPU
            else:
                break
        else:
            if len(null_ibd)<rep:
                print("timed out, killing all processes")
            for p in procs:
                p.terminate()
                p.join()
        for p in procs:
            p.join()
        #except:
            #pass
        while not inqueue.empty():
            inqueue.get()
        for k in cond_null_ibd.keys():
            self.cond_null_ibd[k]=cond_null_ibd[k]
        for x in set(null_ibd):
            null_nd[x]=null_ibd.count(x)/len(null_ibd)
        null_ibd=[]
        all_ibd={}
        return null_nd

    def nulldist_trait(self,n_jobs,rep=5000):
        print "calculating null dist for a single fam.."
        nfnum = len(self.family.nonfounder)
        ###parallel processing###
        inqueue = multiprocessing.Queue()
        manager=Manager()
        null_ibd = manager.list([])
        traits=[self.family.fam_dict[iid]['trait'] for iid in self.family.all_members]
        raw_total=math.factorial(len(self.family.all_members))
        if raw_total > rep:
            traits_full_permutations=[]
            for x in range(rep):
                random.shuffle(traits)
                traits_full_permutations.append(tuple(traits))
        else:
            traits_full_permutations=list(permutations(traits,len(self.family.all_members)))
        for i in traits_full_permutations:
            inqueue.put(i)
        full_size=len(traits_full_permutations)
        #start n seperate processes
        #shuffle phenotypes
        procs = []
        for proc in range(n_jobs):
            p = myProcessTrait(proc,self,inqueue,null_ibd)
            p.start()
            procs.append(p)
            inqueue.put(None)
        TIMEOUT = 3600
        start = time.time()
        while time.time() - start <= TIMEOUT:
            if any([p.is_alive() for p in procs]) and len(null_ibd)<full_size:
                time.sleep(.1)  # Just to avoid hogging the CPU
            else:
                break
        else:
            if len(null_ibd)<full_size:
                print("timed out, killing all processes")
            for p in procs:
                p.terminate()
                p.join()
        for p in procs:
            p.join()
        #except:
        #    pass
        while not inqueue.empty():
            inqueue.get()
        return null_ibd[:]

    def execute(self,m):
        #set the value for matrixs
        self.setmatrix()
        self.combine_inherit_vec(m) #calculate the overall posterior inv probability and corresponding IBD #time-consuming
        self.setcov()      #set covariance matrix #time-consuming
        self.reg()    #regression

def null_generator_pairibd(args):
    data,qnpl=args
    #random.seed(99)
    fgt_rv=[]
    nfnum = len(qnpl.family.nonfounder)
    #calculate null ibd by genotypes
    v,fgt=data
    tmp = bin(v)[2:]                     #get binary format
    inv = '0'*(2*nfnum-len(tmp))+tmp
    for idx,f in enumerate(qnpl.family.founder):
        fgt_rv.extend(['{}{}'.format(f,a_idx) for a_idx,tmp_a in enumerate(fgt[2*idx:2*idx+2]) if tmp_a!=1])
    ibdall=qnpl.cal_ibd_inv(inv,fgt_rv)
    return tuple(ibdall)

def null_generator(procID,founderid,fam,queue,all_ibd,cond_null_ibd,null_ibd,n_jobs,perfect_max):
    #generate permutations by shuffling founder genotypes
    while True:
        try:
            data = queue.get()
            if data is None:
                break
            else:
                foundergt = {}
                new_alleles = {}
                ibdall = [] #pairwise IBD
                nfnum = len(fam.nonfounder)
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
                new_alleles = fam.getGT(inv,foundergt)
                conditional_prob = {}
                #Assign missing again and re-impute
                if fam.missing_all != []:
                    for person in fam.missing_all:#fam.sorted_miss_persons:
                        new_alleles[person]=[0,0]
                allele_key=''
                all_gt=[]
                for person in sorted(new_alleles.keys()):
                    allele_key+=''.join(map(str,new_alleles[person]))
                    all_gt+=new_alleles[person]
                if allele_key not in all_ibd.keys():
                    all_ibd[allele_key]=None
                    if fam.missing_all != [] and set(all_gt)!=set([0,1]):
                        fam.GT_infer(new_alleles,conditional_prob)
                        if fam.rvibd:
                            cond_key=conditional_prob['~combined'][0].items()[0]
                            if cond_key not in cond_null_ibd.keys():
                                cond_null_ibd[cond_key]=None
                                fam.null_permute(n_jobs,perfect_max,conditional_prob)
                                cond_null_ibd[cond_key]=[fam.expect_pair_ibd,fam.prior]
                            else:
                                while cond_null_ibd[cond_key] is None:
                                    time.sleep(.01)
                                fam.expect_pair_ibd,fam.prior=cond_null_ibd[cond_key]
                    #print allele_key
                    #print new_alleles
                    #print all_gt
                    if set(all_gt)==set([0,1]):
                        n,d=None,None
                    else:
                        n,d=fam.matrix_cal(new_alleles,conditional_prob)
                    all_ibd[allele_key]=(n,d)
                    if (n,d)!=(None,None):
                        null_ibd.append((n,d))
                    #print n,d
                else:
                    while all_ibd[allele_key] is None:
                        time.sleep(.01)
                    else:
                        n,d=all_ibd[allele_key]
                        if (n,d)!=(None,None):
                            null_ibd.append((n,d))
        except Exception as e:
            screen_output.err_out("error in null_generator {}:{}".format(procID,e))
            raise TypeError
            break

def null_generator_trait(procID,fam,queue,null_ibd):
    #generate permutations by shuffling phenotypes
    while True:
        try:
            traits = queue.get()
            if traits is None:
                break
            else:
                #traits is one possible permutation of traits
                #permutations
                new_traits={}
                for idx,iid in enumerate(fam.all_members):
                    new_traits[iid]=traits[idx]
                n,d=fam.matrix_cal_trait(new_traits)
                null_ibd.append((n,d))
        except Exception as e:
            screen_output.err_out("error in null_generator {}:{}".format(procID,e))
            raise TypeError
            break

class myProcess(Process):
    def __init__(self, procID, founderid, fam, q, all_ibd, cond_null_ibd,null_ibd,n_jobs,perfect_max):
        Process.__init__(self)
        self.founderid = founderid
        self.q = q
        self.fam = fam
        self.procID = procID
        self.null_ibd = null_ibd
        self.all_ibd = all_ibd
        self.cond_null_ibd = cond_null_ibd
        self.n_jobs = n_jobs
        self.perfect_max = perfect_max
    def run(self):
        try:
            null_generator(self.procID,self.founderid, self.fam, self.q, self.all_ibd, self.cond_null_ibd, self.null_ibd)
        except Exception as e:
            screen_output.err_out("error in: %d"%self.procID)

class myProcessTrait(Process):
    def __init__(self, procID, fam, q, null_ibd):
        Process.__init__(self)
        self.q = q
        self.fam = fam
        self.procID = procID
        self.null_ibd = null_ibd
    def run(self):
        try:
            null_generator_trait(self.procID, self.fam, self.q, self.null_ibd)
        except Exception as e:
            screen_output.err_out("error in: %d"%self.procID)

