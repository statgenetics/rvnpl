#!/usr/bin/python 
#Author: Linhai Zhao
#regression-based QTL mapping in extended families (Sham 2002)
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
from sympy import Symbol,Interval
from RVNPLcpp import cmissingparents,cmissing_infer,cpostInv
from RVNPL import screen_output,utils,family_class_combined,qnpl_class


def execute(args):
        screen_output.log_out(args.path)
        famstruct = []
        famstruct_null = {}
        for (dirpath,dirname,filename) in os.walk(args.path):
                fped = [f for f in filename if re.search(r'ped',f)]
                fdat = [d for d in filename if re.search(r'dat',d)]
                if len(fped) == 1 and len(fdat) == 1:
                    family,mafs,markername,fam_ids=utils.load_data(fped,fdat,dirpath,args,True)
                    tmp_dir=os.path.dirname(dirpath)
                    pdirname = os.path.dirname(tmp_dir)
                    genename=os.path.basename(pdirname)
                    cdist=[]    #combined distribution for each replicate
                    #output original result
                    result = file("{}/original_result.txt".format(args.output),'a')
                    pv = file("{}/pvalue.txt".format(args.output),'a')
                    rtext = '%s:\n'%genename
                    marker_n=len(markername)
                    markers_to_analyze=range(marker_n)
                    numerator = {}                          #marker:[nom1,nom2...]               
                    denominator = {}
                    asymp_pv = {}
                    sum_null_nd = {}
                    Q_all = {}
                    fam_conditional_prob=[]
                    fam_to_analyze={}
                    wt_count={}
                    fam_class={}
                    fam_struct_info = []
                    fam_struct_info_dic = {}
                    fam_struct_cond_null = []
                    fam_struct_cond_nulldic = {}
                    for mkid in markers_to_analyze:
                        fam_to_analyze[markername[mkid]]=copy.copy(fam_ids)
                        fam_class[mkid]={}
                        wt_count[markername[mkid]]=0
                    for m in markers_to_analyze:
                        print markername[m]
                        sum_null_nd[markername[m]]=[]
                        if None in fam_to_analyze[markername[m]]:
                            continue
                        for fid in fam_to_analyze[markername[m]]:
                            print "processing family {} ...".format(fid)
                            fam = family_class_combined.QFamily(m)
                            fam.rvibd = args.rvibd
                            if args.snv and mafs != {}:
                                fam.mp_freq = mafs[markername[m]]
                            elif not args.snv:
                                if mafs[markername[m]].has_key(fid):
                                    fam.mp_freq=mafs[markername[m]][fid]
                            fam.setdict(family[fid])  #set family data structure as a dictionary
                            if fam.wt_fam:
                                wt_count[markername[m]]+=1
                            fam.classify_affect(m) #classify affected individuals into relative pairs and calculate IBD
                            if fam.err == 1 or not fam.info:
                                fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
                                fam.clean()
                                continue
                            else:
                                #get the prior IBD product for each pair of affected relative pairs
                                #fam.prior=sum(prob*pi(i,j)*pi(k,l))
                                fam_qnpl=qnpl_class.QNPL(fam)
                                fam_class[m][fid]=(fam,fam_qnpl)
                                if not args.rvibd:
                                    pfamstruct=copy.copy(fam.famstruct)
                                    pfamstruct.pop('info',None)
                                    if pfamstruct not in famstruct:
                                        #print "new structure"
                                        famstruct.append(pfamstruct)
                                        fam_qnpl.null_inv()
                                        famstruct_null[famstruct.index(pfamstruct)]=fam.prior
                                    else:
                                        idx = famstruct.index(pfamstruct)
                                        fam.prior = famstruct_null[idx]
                                all_member_gt=[]
                                all_member_trait=[]
                                for ind in fam.all_members:
                                    all_member_gt.append(fam.fam_dict[ind]['gt'][2*m:2*m+2])
                                    all_member_trait.append(fam.fam_dict[ind]['trait'])
                                #get the null distribution of numerator and denominator via permutations
                                if fam.missing_all!=[]:
                                    fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                    fam_combined_info_gt=(fam.famstruct,fam.conditional_prob['~combined'][0],all_member_gt)
                                else:
                                    fam_combined_info=(fam.famstruct)
                                    fam_combined_info_gt=(fam.famstruct,all_member_gt)
                                print fam.famstruct
                                #if fam.wt_fam:
                                    #continue
                                #print fam_combined_info
                                if fam_combined_info not in fam_struct_info:
                                    fam_struct_info.append(fam_combined_info)
                                    if args.rvibd:
                                        fam_qnpl.null_permute(n_jobs=args.n_jobs,perfect_max=args.perfect_max)
                                    fam_qnpl.execute(m)
                                    fam_struct_info_dic[fam_struct_info.index(fam_combined_info)]=[fam_qnpl.ES,fam_qnpl.ED,fam_qnpl.ET,fam_qnpl.covSS,fam_qnpl.covDD,fam_qnpl.covSD,fam_qnpl.covY_inverse,fam.prior,fam.expect_pair_ibd]
                                else:
                                    fam_qnpl.ES,fam_qnpl.ED,fam_qnpl.ET,fam_qnpl.covSS,fam_qnpl.covDD,fam_qnpl.covSD,fam_qnpl.covY_inverse,fam.prior,fam.expect_pair_ibd=fam_struct_info_dic[fam_struct_info.index(fam_combined_info)]
                                    fam_qnpl.execute(m)
                                #print fam_combined_info_gt
                                print fam_qnpl.numerator, fam_qnpl.denominator
                                if m not in numerator:
                                    numerator[m] = []
                                    denominator[m] = []
                                numerator[m].append(fam_qnpl.numerator)
                                denominator[m].append(fam_qnpl.denominator)
                                #fam.clean()
                        try:
                            Q = sum(numerator[m])/sum(denominator[m])
                        except:
                            Q = 0
                        t=0
                        if sum(denominator[m])>0:
                            if Q>1:
                                t=sum(denominator[m])
                            elif Q>0:
                                t = Q**2*sum(denominator[m])
                        else:
                            print "Var:",sum(denominator)
                        p_asymp = stats.chi2.sf(t,1)/2
                        Q_all[markername[m]]=(Q,t)
                        asymp_pv[markername[m]]=p_asymp
                    p_emp_min=9
                    apv_min=min(asymp_pv.values())
                    for ele in sorted(asymp_pv.items(), key=lambda x:x[1])[:3]:
                        apv=ele[1]
                        mname=ele[0]
                        idx_apv=markername.index(mname)
                        Q,t=Q_all[mname]
                        p_emp=9
                        if args.exact:
                            fam_reps=[100, 500, 1000, args.fam_rep]
                            reps=[500, 5000, 50000, args.rep]
                            fam_rep=args.fam_rep
                            rep=args.rep
                            if not args.force:
                                if apv>0.5:
                                    fam_rep=fam_reps[0]
                                    rep=reps[0]
                                elif apv>0.1:
                                    fam_rep=fam_reps[1]
                                    rep=reps[1]
                                elif apv>0.01:
                                    fam_rep=fam_reps[2]
                                    rep=reps[2]
                            fam_count=len([x for x in fam_to_analyze[mname] if x is not None])
                            if apv>=args.cut or apv<args.lower_cut:
                                #Adaptive permutations
                                #Use asymptotic p-values if the results are not pointwise significant
                                p_emp=apv
                            while p_emp==9 or p_emp==1/(1+rep):
                                #permutations needed for empirical p-values
                                #increase_permutation=False
                                if p_emp==1/(1+rep):
                                    #if empiric p-value reached limit, more permutation is needed
                                    if fam_rep==args.fam_rep:
                                        break
                                    #increase_permutation=True
                                    sum_null_nd[mname]=[]
                                    print "increase permutation"
                                    print fam_rep,rep
                                    idx=fam_reps.index(fam_rep)
                                    fam_rep=fam_reps[idx+1]
                                    rep=reps[idx+1]
                                print fam_rep,rep
                                for fid in fam_to_analyze[mname]:
                                    if fid != None:
                                        print fid
                                        fam,fam_qnpl = fam_class[idx_apv][fid]
                                        if args.pheno:
                                            sum_null_nd[mname].append(fam_qnpl.nulldist_trait(n_jobs=args.n_jobs,rep=fam_rep))
                                        else:
                                            if fam.wt_fam:
                                                null_nd={}
                                                null_nd[(fam_qnpl.numerator,fam_qnpl.denominator)]=1
                                            else:
                                                if fam.missing_all!=[]:
                                                    fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                                    fam_combined_info_gt=(fam.famstruct,fam.conditional_prob['~combined'][0],all_member_gt)
                                                else:
                                                    fam_combined_info=(fam.famstruct)
                                                    fam_combined_info_gt=(fam.famstruct,all_member_gt)
                                                if fam_combined_info not in fam_struct_cond_null:
                                                    fam_struct_cond_null.append(fam_combined_info)
                                                    null_nd=fam_qnpl.nulldist(n_jobs=args.n_jobs,perfect_max=args.perfect_max,rep=fam_rep)
                                                    fam_struct_cond_nulldic[fam_struct_cond_null.index(fam_combined_info)]=fam_qnpl.cond_null_ibd
                                                else:
                                                    fam_qnpl.cond_null_ibd=fam_struct_cond_nulldic[fam_struct_cond_null.index(fam_combined_info)]
                                                    null_nd=fam_qnpl.nulldist(n_jobs=args.n_jobs,perfect_max=args.perfect_max,rep=fam_rep)
                                            sum_null_nd[mname].append(null_nd)
                                cdist = utils.randomsample_qtl(sum_null_nd[mname],rep=rep)
                                print cdist[:10], cdist[-10:]
                                print Q,t
                                p_emp=(len([q for q in cdist if q >= t])+1-cdist.count(t)/2)/(len(cdist)+1)
                        rtext += '%s\t%.4f\t%.4f\t%.7f\t%.7f\n'%(mname,Q,t,apv,p_emp)
                        p_emp_min = p_emp if p_emp<p_emp_min else p_emp_min
                    result.write(rtext)
                    if args.exact:
                        pv.write('%s\t%.7f\n'%(genename,p_emp_min))
                    else:
                        pv.write('%s\t%.7f\n'%(genename,apv_min))
                    pv.close()
                    result.close()
