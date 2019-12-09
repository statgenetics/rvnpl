#!/usr/bin/python
from __future__ import division
from multiprocessing import Process, Value, Manager, Pool
from multiprocessing.sharedctypes import Array
import multiprocessing
import random
import copy
from ctypes import c_double
from itertools import product
import math
import numpy as np
import time
from scipy import stats
from scipy.optimize import fsolve
from sympy.solvers import solve,solveset
from sympy import Symbol,Interval
from RVNPLcpp import sall_cpp,cpostInv
from RVNPL import screen_output,ibd_module
manager=Manager()
class NPL:
        def __init__(self,family):
            self.family=family
            self.ibd_total=0
            self.ibd_sall=0 
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

        def sall(self,v,fam_allele,fmarker):
            ##calculate S-all based on given inheritance vector#####
            h = []           #collections of alleles from affected individuals
            try:
                self.family.affnf[0]+1
                s_all = sall_cpp.apply(self.family.affected,fam_allele,self.family.fallele,fmarker,self.family.rvibd)
            except TypeError:
                #code affected if not int
                code_dic={}
                code_int=0
                for tmp_af in self.family.affected:
                    code_int+=1
                    code_dic[tmp_af]=code_int
                tmp_aff=[code_dic[tmp_af] for tmp_af in self.family.affected]
                tmp_fam_allele={}
                for tmp_key in self.family.affected:
                    tmp_fam_allele[code_dic[tmp_key]]=fam_allele[tmp_key]
                if fmarker==[]:
                    s_all = sall_cpp.apply(tmp_aff,tmp_fam_allele,self.family.fallele,fmarker,False)
                else:
                    s_all = sall_cpp.apply(tmp_aff,tmp_fam_allele,self.family.fallele,fmarker,self.family.rvibd)
            return s_all

        def cal_sall(self,pinv,fmarker,pinv_sall=[],pinv_key_dict={},tmp_alleles={},t_alleles={},rvibd=False):
            tmp_sall=0
            ###NPL-all
            if rvibd: #if RV only, sall is dependent on founder alleles and post inv
                aff_rv_alleles=[]  #inheritance code for RVs in affecteds
                dup_person={}
                multi_a={}
                for af in self.family.affected:
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
                    for afid in self.family.affected:
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
                    for afid in self.family.affected:
                            aff_allele_pool.extend(t_alleles[afid])
                    for founder_a in list(set(aff_allele_pool)):
                        fa_count=aff_allele_pool.count(founder_a)
                        if fa_count>1:
                            tmp_code.append(fa_count)
                            multi_a[founder_a]=[fa_count,tmp_code.count(fa_count)]
                    for afid in self.family.affected:
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
                    conditional_prob=self.family.conditional_prob
                if tmp_alleles=={}:
                    for iid in self.family.fam_dict.keys():
                        tmp_alleles[iid]=self.family.fam_dict[iid]['gt'][2*self.family.mid:2*self.family.mid+2]
                if pinv_pair_dict == {}:
                    pinv_pair_dict = self.pinv_pair_dict
                if pinv_sall == []:
                    pinv_sall = self.pinv_sall
                    if self.pinv_sall == []:
                        self.pinv_sall = Array(c_double,[0 for x in xrange(pow(2,2*len(self.family.invnf)))])
                        pinv_sall=self.pinv_sall
                if pinv_key_dict == {}:
                    pinv_key_dict = self.pinv_key_dict
                tmp_dic={}
                gt_id=[]
                output_ibd=0
                output_sall=0
                if self.family.sorted_miss_persons != []:
                    tmp_dic=conditional_prob['~combined'][0]
                    gt_id=conditional_prob['~combined'][1]
                else:
                    tmp_dic={'0000':1}  #dummy dic when no missing
                for tmp_gt,prob in tmp_dic.iteritems():
                    total_ibd = 0
                    total_sall = 0
                    #print tmp_gt,prob
                    #fill out missing genotypes based on inference
                    for m_person in self.family.sorted_miss_persons:
                        idx=gt_id.index(m_person[0])
                        tmp_alleles[m_person[0]]=[ord(x)-96 for x in tmp_gt[2*idx:2*idx+2]]
                    if self.family.rvibd:
                        cal_flag=False
                        affmarker=[]
                        if sall_flag:
                            for afid in self.family.affected:
                                affmarker+=tmp_alleles[afid]
                            for uniq_marker in list(set(affmarker)):
                                if uniq_marker!=1 and affmarker.count(uniq_marker)>1:
                                    cal_flag=True
                                    break
                        else:
                            aff_ind=list(set([tmp_aff for pair in self.family.affected_pairs for tmp_aff in pair]))
                            for afid in aff_ind:
                                affmarker+=tmp_alleles[afid]
                            for uniq_marker in list(set(affmarker)):
                                if uniq_marker!=1 and affmarker.count(uniq_marker)>1:
                                    tmp_share_pattern={}
                                    for afid in aff_ind:
                                        tmp_share_pattern[afid]=[True,True] if uniq_marker in tmp_alleles[afid] else [False,False]
                                        if self.family.generation(afid)<-1 and afid not in self.family.founder and uniq_marker in tmp_alleles[afid]:
                                            nf_parent=[tmp_parent for tmp_parent in self.family.parents[afid] if tmp_parent not in self.family.founder][0]
                                            f_parent=self.family.parents[afid][1-self.family.parents[afid].index(nf_parent)]
                                            if uniq_marker not in tmp_alleles[nf_parent] and f_parent not in self.family.affected:
                                                tmp_share_pattern[afid][1]=False
                                    for tmp_pair in self.family.affected_pairs:
                                        if tmp_share_pattern[tmp_pair[0]][0] and tmp_share_pattern[tmp_pair[1]][0]:
                                            if self.family.parents[tmp_pair[0]]==self.family.parents[tmp_pair[1]]:
                                                cal_flag=True
                                                break
                                            elif tmp_share_pattern[tmp_pair[0]][1] or tmp_share_pattern[tmp_pair[1]][1]:
                                                cal_flag=True
                                                break
                                    if cal_flag:
                                        break
                        if not cal_flag:
                            total_ibd = 0
                            if sall_flag:
                                total_sall = 1
                                output_sall+=total_sall*prob
                            continue
                    fmarker=[]
                    for f in self.family.invfounder:
                        fmarker+=tmp_alleles[f]
                    if self.family.simple:
                        #NPL-pair
                        #print self.sib, self.cousin, self.un
                        gpid,gmid=self.family.gpid,self.family.gmid
                        ibd_n=0        #number of alleles IBD among affected individuals        
                        for aid,bid in self.family.offspring_pairs:
                            if not self.family.rvibd:
                                total_ibd+=1
                            else:
                                if tmp_alleles[aid]==[1,1] or tmp_alleles[bid]==[1,1]:
                                    continue
                                else:
                                    other_parent=[tmp_p for tmp_p in self.family.fam_dict[aid]['parents'] if tmp_p != bid][0]
                                    genotype = tmp_alleles[aid]+tmp_alleles[bid]+tmp_alleles[other_parent]
                                    total_ibd+=ibd_module.child_ibd(genotype)
                        for aid,bid in self.family.sib:
                            ap = self.family.fam_dict[aid]['parents']
                            genotype = tmp_alleles[ap[0]]+tmp_alleles[ap[1]]+tmp_alleles[aid]+tmp_alleles[bid]
                            ibd = ibd_module.sib_ibd(genotype,rv_flag=self.family.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.family.cousin:
                            ap = copy.copy(self.family.fam_dict[aid]['parents'])
                            bp = copy.copy(self.family.fam_dict[bid]['parents'])
                            self.family.firstparent(ap,gpid,gmid)
                            self.family.firstparent(bp,gpid,gmid)
                            genotype = tmp_alleles[gpid]+tmp_alleles[gmid]+tmp_alleles[ap[0]]+\
                                    tmp_alleles[ap[1]]+tmp_alleles[aid]+tmp_alleles[bp[0]]+\
                                    tmp_alleles[bp[1]]+tmp_alleles[bid]
                            ibd = ibd_module.cousin_ibd(genotype,rv_flag=self.family.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.family.hsib:
                            ap = copy.copy(self.family.fam_dict[aid]['parents'])
                            bp = copy.copy(self.family.fam_dict[bid]['parents'])
                            sp = list(set(ap)&set(bp))[0]                #shared parent between half-sibs
                            ap.remove(sp)
                            bp.remove(sp)
                            genotype = tmp_alleles[sp]+tmp_alleles[ap[0]]+tmp_alleles[aid]+\
                            tmp_alleles[bp[0]]+tmp_alleles[bid]               #[shared parent, non-share parent, kid]
                            ibd=ibd_module.hsib_ibd(genotype,rv_flag=self.family.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.family.un:      #(uncle,nephew)
                            ap = copy.copy(self.family.fam_dict[aid]['parents'])
                            bp = copy.copy(self.family.fam_dict[bid]['parents'])
                            self.family.firstparent(bp,gpid,gmid)
                            genotype = tmp_alleles[gpid]+tmp_alleles[gmid]+tmp_alleles[aid]+\
                                    tmp_alleles[bp[0]]+tmp_alleles[bp[1]]+tmp_alleles[bid]
                            ibd = ibd_module.un_ibd(genotype,rv_flag=self.family.rvibd)
                            total_ibd+=ibd
                        for aid,bid in self.family.gp:      #(grandpa,grandson)
                            bp = copy.copy(self.family.fam_dict[bid]['parents'])
                            self.family.firstparent(bp,gpid,gmid)
                            grand_other=gmid if aid==gpid else gpid
                            genotype = tmp_alleles[aid]+tmp_alleles[grand_other]+\
                                    tmp_alleles[bp[0]]+tmp_alleles[bp[1]]+tmp_alleles[bid]
                            ibd = ibd_module.gp_ibd(genotype,rv_flag=self.family.rvibd)
                            total_ibd+=ibd
                        if sall_flag:
                            ###NPL-all
                            #get posterior inheritance vectors for probable allele configurations
                            postinv=cpostInv.apply(self.family.invnf, self.family.founder, self.family.parents, tmp_alleles, self.family.rvibd)
                            len_postinv=len(postinv)
                            #print fmarker
                            #print "postinv:"+repr(len_postinv)+" pinv_key_dict:"+repr(len(pinv_key_dict.keys()))
                            postinv_list=set(postinv) if self.family.rvibd else postinv
                            for pinv in postinv_list:
                                t_alleles=self.family.assign_allele(pinv)
                                t_sall=self.cal_sall(pinv,fmarker,pinv_sall,pinv_key_dict,tmp_alleles,t_alleles,rvibd=self.family.rvibd)
                                if self.family.rvibd:
                                    total_sall+=t_sall*postinv.count(pinv)
                                else:
                                    total_sall+=t_sall
                            total_sall=total_sall/len_postinv
                    else:
                        if self.family.rvibd:
                            founder_rv_alleles=[]
                            for founder in self.family.founder:
                                for tmp_i in range(2):
                                    if tmp_alleles[founder][tmp_i] not in [0,1]:
                                        #alt allele, i.e. RV
                                        founder_rv_alleles.append('{}{}'.format(founder,tmp_i))
                        #get posterior inheritance vectors for probable allele configurations
                        postinv=cpostInv.apply(self.family.invnf, self.family.founder, self.family.parents, tmp_alleles, self.family.rvibd)
                        len_postinv=len(postinv)
                        #print "postinv:"+repr(len_postinv)+" pinv_pair:"+repr(len(pinv_pair_dict.keys()))
                        postinv_list=set(postinv) if self.family.rvibd else postinv
                        for pinv in postinv_list:
                            t_alleles={}
                            ###NPL-pairs
                            int_pinv=int(pinv,2)
                            shared_alleles=[]
                            t_ibd=0
                            if int_pinv in pinv_pair_dict:
                                shared_alleles=pinv_pair_dict[int_pinv]
                            else:
                                t_alleles=self.family.assign_allele(pinv)
                                for af_pair in self.family.affected_pairs:
                                    for tt_a in t_alleles[af_pair[0]]:
                                        if tt_a in t_alleles[af_pair[1]]:
                                            shared_alleles.append(tt_a)
                                    #tmp_ibd=set(t_alleles[af_pair[0]])&set(t_alleles[af_pair[1]])
                                    #shared_alleles.extend(list(tmp_ibd))
                                pinv_pair_dict[int_pinv]=shared_alleles
                            if self.family.rvibd:
                                t_ibd=len([tmp_sh_a for tmp_sh_a in shared_alleles if tmp_sh_a in founder_rv_alleles])
                                total_ibd+=t_ibd*postinv.count(pinv)
                            else:
                                t_ibd=len(shared_alleles)
                                total_ibd+=t_ibd
                            if sall_flag:
                                ###NPL-all
                                if t_alleles=={}:
                                    t_alleles=self.family.assign_allele(pinv)
                                t_sall=self.cal_sall(pinv,fmarker,pinv_sall,pinv_key_dict,tmp_alleles,t_alleles,self.family.rvibd)
                                if self.family.rvibd:
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

        def nullibd(self,rep,n_jobs,sall_flag=False,infer_flag=2,simple=False,verbose=1):
            #calculate expected mean and std for IBD under H0 
            if verbose==1:
                screen_output.run_out("calculating null distribution...")
            #the number of nonfounders that should be included in inheritance vector
            founderid = self.family.founder
            ###parallel processing###
            inqueue = multiprocessing.Queue()
            all_ibd = manager.dict(self.all_ibd)
            null_ibd = manager.list(self.null_ibd)
            if self.pinv_sall==[]:
                pinv_sall = Array(c_double,[0 for x in xrange(pow(2,2*len(self.family.invnf)))])
            else:
                pinv_sall=(c_double*len(self.pinv_sall))(*self.pinv_sall)
            pinv_key_dict = manager.dict(self.pinv_key_dict)
            pinv_pair_dict = manager.dict(self.pinv_pair_dict)
            for i in xrange(rep):
                inqueue.put(i)
            #start n seperate processes
            procs = []
            for proc in range(n_jobs):
                p = myProcess(proc,self.family,self,inqueue,all_ibd,\
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

        def null_perfect(self,n_jobs,perfect_max,update_null=True,verbose=1):
        #null distribution of IBDs under perfect data
            if verbose==1:
                screen_output.run_out("calculating null distribution using perfect data approximation...")
            p = Pool(processes=n_jobs)
            inv_list=[]
            rep=pow(2,2*len(self.family.invnf))
            #print rep
            if rep>perfect_max:
                inv_list=random.sample(xrange(rep),perfect_max)
            else:
                inv_list=range(rep)
            self.pinv_sall = []
            chunk_size=len(inv_list)//n_jobs if len(inv_list)>=n_jobs else 1
            result = p.imap(null_inv, product(inv_list,[(self.family,self)]), chunksize=chunk_size)
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

        def null_perfect_rvibd(self,n_jobs,perfect_max,sall_flag=False,infer_flag=2,verbose=1):
        #null distribution of IBDs under perfect data
            if verbose==1:
                screen_output.run_out("calculating theoretical null distribution for rvibd...")
            rep=pow(2,2*len(self.family.nonfounder))
            #print rep
            #Get all possible parental genotype configurations
            foundergt_list=[]
            foundergt_dic={}
            for f in self.family.founder:
                foundergt_list += self.family.fam_dict[f]['gt'][2*self.family.mid:2*self.family.mid+2]
            try:
                #if there are missing
                gt_id=self.family.conditional_prob['~combined'][1]
                for gt,prob in self.family.conditional_prob['~combined'][0].iteritems():
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
            mean_pair, mean_all, var_pair, var_all = 0, 0, 0, 0
            all_ibd = manager.dict(self.all_ibd)
            if self.pinv_sall==[]:
                pinv_sall = Array(c_double,[0 for x in xrange(pow(2,2*len(self.family.invnf)))])
            else:
                pinv_sall=(c_double*len(self.pinv_sall))(*self.pinv_sall)
            pinv_key_dict = manager.dict(self.pinv_key_dict)
            pinv_pair_dict = manager.dict(self.pinv_pair_dict)
            combined_dist={}
            prob_and_mean=[]
            if self.family.simple:
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
                    p = myProcess(proc,self.family,self,inqueue,all_ibd,\
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

def null_inv(args):
#calculate null ibd by given inheritance vector
    v,(fam,fam_npl)=args
    tmp = bin(v)[2:]                     #get binary format
    inv = '0'*(2*len(fam.invnf)-len(tmp))+tmp
    inv = [int(tmpv) for tmpv in inv]
    fam_alleles=fam.assign_allele(inv)
    ibd_n=0 #sum of pairwise IBD among affected pairs
    for aff_pair in fam.affected_pairs:
        ibd_pair=len(set(fam_alleles[aff_pair[0]])&set(fam_alleles[aff_pair[1]]))
        ibd_n+=ibd_pair
    tmp_sall=fam_npl.sall(inv,fam_alleles,[])
    return ibd_n,tmp_sall,v


def null_generator(procID,fam,fam_npl,queue,all_ibd,null_ibd,pinv_sall,pinv_key_dict,pinv_pair_dict,sall_flag,infer_flag):
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
                    for f in fam.founder:
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
                    ibd_result=fam_npl.cal_ibd(new_alleles,conditional_prob,pinv_pair_dict,pinv_sall,pinv_key_dict,sall_flag)
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
    def __init__(self, procID, fam, fam_npl, q,all_ibd,null_ibd,\
            pinv_sall,pinv_key_dict,pinv_pair_dict,sall_flag,infer_flag=2):
        Process.__init__(self)
        self.q = q
        self.fam = fam
        self.fam_npl = fam_npl
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
            null_generator(self.procID, self.fam, self.fam_npl,\
                    self.q, self.all_ibd,self.null_ibd,\
                    self.pinv_sall,self.pinv_key_dict,self.pinv_pair_dict,self.sall_flag,self.infer_flag)
        except Exception as e:
            screen_output.err_out("error in: %d"%self.procID)
