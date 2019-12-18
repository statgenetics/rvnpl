#!/usr/bin/python 
#Author: Linhai Zhao
#NPL scoring functions to detect over-sharing in affecteds in extended families
from __future__ import division
#from memory_profiler import profile
import os
import re
import random
import copy
from zipfile import ZipFile
from ctypes import *
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
from RVNPL import screen_output,utils,family_class_combined,npl_class
#################
#manager=Manager()

#############################################
def execute(args):
        if args.verbose>=0:
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
                    family,mafs,markername,fam_ids=utils.load_data(fped,fdat,dirpath,args)
		    tmp_dir=os.path.dirname(dirpath)
		    pdirname = os.path.dirname(tmp_dir)
		    genename=os.path.basename(pdirname)#foldername.split('_gnomAD')[0]
		    #output original result
		    if args.output is None:
			args.output = pdirname
		    result = file("{}/original_result.txt".format(args.output),'a')
		    rtext = '%s:\n'%genename
		    #take in ped files
		    fam_class = {}
		    marker_n=len(markername)
		    fam_num = len(family)
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
                        if args.verbose >=0:
			    screen_output.run_out("analyzing "+repr(markername[m]))
			if None in fam_to_analyze[markername[m]]:
			    continue
			for fid in fam_to_analyze[markername[m]]:
			    #for each family, construct a class of Family
                            if args.verbose==1:
			        print "fid:"+repr(fid)
			    fam = family_class_combined.Family(m)
			    fam.rvibd = args.rvibd
			    if args.snv and mafs != {}:
				fam.mp_freq = mafs[markername[m]]
			    elif not args.snv:
				if mafs[markername[m]].has_key(fid):
				    fam.mp_freq=mafs[markername[m]][fid]
			    fam.setdict(family[fid])                   #set family dictionary
			    try:
				if args.info_only and not fam.info:
                                    if args.verbose==1:
				        screen_output.run_out("uninformative family detected")
				    fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
				    fam.clean()
				    continue
				fam.classify_affect()
			    except:
				fam.err=True
			    if fam.err:
                                if args.verbose==1:
				    screen_output.err_out("Error detected in fam {} for marker {}".format(fid,markername[m]))
				fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
				fam.clean()
				continue
                            fam_npl=npl_class.NPL(fam)
			    fam_class[m][fid]=(fam,fam_npl)
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
                                            fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=args.sall,verbose=args.verbose)
                                            fam_conditional_prob_perfect.append(fam_combined_info)
                                            fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                        else:
                                            fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std=fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]
                                    else:
                                        if fam.famstruct not in famstruct_perfect:
                                            fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=args.sall,verbose=args.verbose)
                                            famstruct_perfect.append(fam.famstruct)
                                            famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                        else:
                                            fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]
                                else:
                                    pfamstruct=copy.copy(fam.famstruct)
                                    pfamstruct.pop('info',None)
                                    if pfamstruct not in famstruct_perfect:
                                        fam_npl.null_perfect(n_jobs=args.n_jobs,perfect_max=args.perfect_max,verbose=args.verbose)
                                        famstruct_perfect.append(pfamstruct)
                                        famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]=[fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                    else:
                                        fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]
			    else:
				#use permutations to calculate nullibd
				if fam.missing_all!=[] and not fam.wt_fam:
                                    fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                    if fam_combined_info not in fam_conditional_prob:
                                        fam_npl.nullibd(rep=50,n_jobs=args.n_jobs,sall_flag=args.sall,verbose=args.verbose)
                                        fam_conditional_prob.append(fam_combined_info)
                                        fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=[fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                    else:
                                        fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std=fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]
				elif fam.info and not fam.wt_fam:   
				    if fam.famstruct not in famstruct:
                                        fam_npl.nullibd(rep=100,n_jobs=args.n_jobs,sall_flag=args.sall,verbose=args.verbose)
                                        famstruct.append(fam.famstruct)
                                        fam_npl.dist_s=fam_npl.distribution(pall_flag=args.sall)
                                        famstruct_null[str(famstruct.index(fam.famstruct))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                    else:
                                        fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std=famstruct_null[str(famstruct.index(fam.famstruct))]
			    ###############
			    allele_key=''
			    for person in sorted(fam.fam_dict.keys()):
				if person in fam.missing_all:
				    allele_key+='00'
				else:
				    allele_key+=''.join(map(str,fam.fam_dict[person]['gt'][2*fam.mid:2*fam.mid+2]))
			    if allele_key in fam_npl.all_ibd:
				fam_npl.ibd_total, fam_npl.ibd_sall=fam_npl.all_ibd[allele_key]
			    else:
				try:
				    fam_npl.ibd_total,fam_npl.ibd_sall=fam_npl.cal_ibd(sall_flag=args.sall)
				except:
				    fam.err=True
			    if fam.err:
                                if args.verbose==1:
				    screen_output.err_out("Error detected in fam {} for marker {}".format(fid,markername[m]))
				fam_to_analyze[markername[m]][fam_ids.index(fid)]=None
				fam.clean()
				continue
			    z,o = fam_npl.stat()                             #get Z-score and observed IBD among affected individuals for each family 
			    if not fam.wt_fam:
				z_sum[m] += z
			    z_pair[m].append(z)
			    observe[m] = float('%.9f'%(observe[m]+o))                     #combined value of IBD across family
                            if args.sall:
                            ###S-all statistic########
				z_sall=0
				try:
				    z_sall=(fam_npl.ibd_sall-fam_npl.sall_null_mean)/fam_npl.sall_null_std
                                except:
				    z_sall=0
				if not fam.wt_fam:
				    zall_sum[m] += z_sall
				z_all[m].append(z_sall)
                                ###############
                            if args.perfect:
                                if args.kc:
                                    b_pair=fam_npl.null_std/(fam_npl.null_mean-min([x[0] for x in fam_npl.null_ibd]))
                                    prior_b_pair[m].append(b_pair)
                                    if args.sall:
                                        min_sall=min([x[1] for x in fam_npl.null_ibd])
                                        b_all=fam_npl.sall_null_std/(fam_npl.sall_null_mean-min_sall)
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
			asymp_pv[m] = (p_norm,p_norm)
			z_scores[m] = (Z_pair,Z_pair)
			if args.sall:
			    #asymptotic p-value NPL-all
			    Z_all = zall_sum[m]/math.sqrt(fam_num)
			    pall_norm = stats.norm.sf(Z_all)
			    asymp_pv[m] = (p_norm,pall_norm)
			    z_scores[m] = (Z_pair,Z_all)
			if args.kc and z_scores[m][0]>0 and z_scores[m][1]>0:
			    ####Kong&Cox#####
			    ##only tested for delta>0###
			    coeff_pair=[]
			    for tmp_z_pair in z_pair[m]:
				coeff_pair.append(tmp_z_pair/math.sqrt(fam_num))
			    min_b_pair=min(prior_b_pair[m])*math.sqrt(fam_num)
			    Zpair_kc,p_pair_norm_kc=utils.K_C(coeff_pair,min_b_pair)	
			    asymp_pv[m] = (p_pair_norm_kc,p_pair_norm_kc)
			    z_scores[m] = (Zpair_kc,Zpair_kc)
			    if args.sall:
				coeff_all=[]
				for tmp_z_all in z_all[m]:
				    coeff_all.append(tmp_z_all/math.sqrt(fam_num))
				min_b_all=min(prior_b_all[m])*math.sqrt(fam_num)
				Zall_kc,p_all_norm_kc=utils.K_C(coeff_all,min_b_all)	
				asymp_pv[m] = (p_pair_norm_kc,p_all_norm_kc)
				z_scores[m] = (Zpair_kc,Zall_kc)
		    pr_s_min, pr_sall_min = 1, 1
		    small_num=3
                    min_asymp_pv = min([x for x in asymp_pv if x !=1])[0]
                    if min_asymp_pv>0.5 or min_asymp_pv<args.lower_cut:
                        small_num=1
		    for tapv in sorted([(x,idx_apv) for idx_apv,x in enumerate(asymp_pv) if x !=1])[:small_num]:
			#take the 3 smallest asymptotic p-value to calculate empiric p-value
			apv=tapv[0]
			idx_apv = tapv[1]
			asymp_z_pair,asymp_z_all=z_scores[idx_apv]
			rtext += 'marker: %s\n'%markername[idx_apv]
			rtext += 'asymptotic p-value:'+repr(apv)+'\n'
			#rtext += 'Z-pair:{}\tZ-all:{}\n'.format(asymp_z_pair,asymp_z_all)
			fam_num = len([x for x in fam_to_analyze[markername[idx_apv]] if x != None])
			pr_s, pr_sall = None, None
                        if args.verbose>=0:
			    print "marker:%s"%markername[idx_apv]
			    print apv
			rep=args.rep
			fam_rep=args.fam_rep
			fam_reps=[50, 100, 1000, args.fam_rep]
                	reps=[500, 5000, 50000, args.rep]
			if fam_num==0:
			    pr_s, pr_sall = 0.5, 0.5
                        if apv[0]>=args.cut and apv[1]>=args.cut or apv[0]<args.lower_cut and apv[1]<args.lower_cut:
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
                            z_sum_origin=0
                            zall_sum_origin=0
			    increase_permutation=False
			    if pr_s==1/(1+rep) or pr_sall==1/(1+rep):
				#if empiric p-value reached limit, more permutation is needed
				if fam_rep==args.fam_rep:
				    break	
				increase_permutation=True
				idx=fam_reps.index(fam_rep)
				fam_rep=fam_reps[idx+1]
				rep=reps[idx+1]
			    pall_flag=True if fam_rep>1000 or args.sall else False
			    infer_full_flag=1 if fam_rep>1000 else 2
			    rtext+='permutations:{}\n'.format(rep)
			    for fid in fam_to_analyze[markername[idx_apv]]:
				if fid != None:
				    fam,fam_npl = fam_class[idx_apv][fid]
				    if not increase_permutation:
					fam_npl.null_ibd = []
				    if not fam.simple:
					infer_full_flag=2
				    if infer_full_flag==1 and fam.missing_all!=[] and not fam.wt_fam:
					fam.GT_infer(infer_flag=infer_full_flag)
					fam_npl.ibd_total,fam_npl.ibd_sall=fam_npl.cal_ibd(sall_flag=args.sall)
					fam_npl.all_ibd={}
				    if fam.missing_all!=[]:
                                        fam_combined_info=(fam.famstruct,fam.conditional_prob['~combined'][0])
                                        if args.perfect and not args.rvibd:
                                            if fam_combined_info not in fam_conditional_prob:
                                                fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True,verbose=args.verbose)
                                                fam_conditional_prob.append(fam_combined_info)
                                                fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=fam_npl.null_ibd
                                            else:
                                                fam_npl.null_ibd=fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]
                                                if pall_flag and fam_npl.null_ibd[-1][1]==0 or len(fam_npl.null_ibd)<fam_rep:
                                                    fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True,verbose=args.verbose)
                                                    fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=fam_npl.null_ibd
					elif args.rvibd:
					    if args.perfect:
						fam_npl.null_ibd=[]
						try:
						    fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std=fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]
						except:
                                            	    fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
						    fam_conditional_prob_perfect.append(fam_combined_info)
						    fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
						if pall_flag and fam_npl.sall_null_mean==0:
                                            	    fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
						    fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
					    elif fam_npl.null_std==0:
						#unlikely to get different values through permutations
						#this can only occur in rvibd when missing
						fam_npl.null_ibd=[]
						if fam_combined_info not in fam_conditional_prob_perfect:
                                            	    fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
						    fam_conditional_prob_perfect.append(fam_combined_info)
						    fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
						else:
						    fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std=fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]
						    if pall_flag and fam_npl.sall_null_mean==0:
                                            	        fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
							fam_conditional_prob_null_perfect[str(fam_conditional_prob_perfect.index(fam_combined_info))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                            else: #use permutations
						fam_npl.dist_s=None
                                                if fam_combined_info not in fam_conditional_prob:
                                                    fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
                                                    fam_conditional_prob.append(fam_combined_info)
                                                    fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=[fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                                else:
                                                    fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std=fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]
                                                    if pall_flag and fam_npl.null_ibd[-1][1]==0 or len(fam_npl.null_ibd)<fam_rep:
                                                        fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
                                                        fam_conditional_prob_null[str(fam_conditional_prob.index(fam_combined_info))]=[fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                        if fam_npl.null_ibd:
                                            add_time=0
                                            while max([x[0] for x in fam_npl.null_ibd])<fam_npl.ibd_total:
                                                if args.perfect and not args.rvibd:
                                                    fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True,verbose=args.verbose)
                                                else:
                                                    fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
                                                add_time+=1
                                                if add_time > 10:
                                                    break
				    elif fam.info:
					if args.perfect and not args.rvibd:
					    if fam.famstruct not in famstruct:
						fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True,verbose=args.verbose)
						famstruct.append(fam.famstruct)
						fam_npl.dist_s=fam_npl.distribution(pall_flag=pall_flag)
                                                famstruct_null[str(famstruct.index(fam.famstruct))]=[fam_npl.dist_s,len(fam_npl.null_ibd)]
					    else:
                                                fam_npl.dist_s,fam_npl.null_ibd_len=famstruct_null[str(famstruct.index(fam.famstruct))]
                                                if pall_flag and fam_npl.dist_s[1][0][1]==0 or fam_npl.null_ibd_len<fam_rep:
                                                    fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,simple=True,verbose=args.verbose)
                                                    fam_npl.dist_s=fam_npl.distribution(pall_flag=pall_flag)
                                                    famstruct_null[str(famstruct.index(fam.famstruct))]=[fam_npl.dist_s,len(fam_npl.null_ibd)]
					else:
                                            if args.perfect:
                                                fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]
						if pall_flag and fam_npl.sall_null_mean==0:
                                            	    fam_npl.null_perfect_rvibd(n_jobs=args.n_jobs,perfect_max=args.perfect_max,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
						    famstruct_null_perfect[str(famstruct_perfect.index(fam.famstruct))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std] 
					    else:
                                                if fam.famstruct not in famstruct:
                                                    fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
                                                    famstruct.append(fam.famstruct)
                                                    fam_npl.dist_s=fam_npl.distribution(pall_flag=pall_flag)
                                                    famstruct_null[str(famstruct.index(fam.famstruct))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]
                                                else:
                                                    fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std=famstruct_null[str(famstruct.index(fam.famstruct))]
                                                    if pall_flag and fam_npl.dist_s[1][0][1]==0 or len(fam_npl.null_ibd)<fam_rep:
                                                        fam_npl.nullibd(rep=fam_rep,n_jobs=args.n_jobs,sall_flag=pall_flag,infer_flag=infer_full_flag,verbose=args.verbose)
                                                        fam_npl.dist_s=fam_npl.distribution(pall_flag=pall_flag)
                                                        famstruct_null[str(famstruct.index(fam.famstruct))]=[fam_npl.dist_s,fam_npl.null_mean,fam_npl.null_std,fam_npl.null_ibd,fam_npl.sall_null_mean,fam_npl.sall_null_std]				
                                    original_dist=[fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
				    if args.rvibd:
					pfamstruct=copy.copy(fam.famstruct)
					pfamstruct.pop('info',None)
					if pfamstruct not in famstruct_perfect:
                                            fam_npl.null_perfect(n_jobs=args.n_jobs,perfect_max=args.perfect_max,update_null=False,verbose=args.verbose)
					    famstruct_perfect.append(pfamstruct)
					    famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]=[fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std]
					else:
					    fam_npl.null_mean,fam_npl.null_std,fam_npl.sall_null_mean,fam_npl.sall_null_std=famstruct_null_perfect[str(famstruct_perfect.index(pfamstruct))]
				    if pall_flag:
					if fam_npl.ibd_sall == 0:
					###S-all statistic########
					    fam_npl.ibd_sall=fam_npl.cal_ibd({},fam.conditional_prob,sall_flag=pall_flag)[1]	
				    tmp = {}
				    if args.exact and cz_dist == []:
					d,v = fam_npl.dist_s if fam_npl.dist_s else fam_npl.distribution(pall_flag=pall_flag)            #get the expected IBD distribution
					if fam_npl.null_mean:
					    try:
						zv_pair = (fam_npl.ibd_total-fam_npl.null_mean)/fam_npl.null_std
					    except ZeroDivisionError:
						zv_pair = 0
					else:
					    zv_pair = 0
					z_sum_precise+=float('%.9f'%zv_pair)
                                        try:
					    z_sum_origin+=float('%.9f'%((fam_npl.ibd_total-original_dist[0])/original_dist[1]))
                                        except:
                                            pass
					if pall_flag:
					    if fam_npl.sall_null_mean:
						try:
						    zv_all = (fam_npl.ibd_sall-fam_npl.sall_null_mean)/fam_npl.sall_null_std
						except ZeroDivisionError:
						    zv_all = 0
					    else:
						zv_all = 0
					    zall_sum_precise+=float('%.9f'%zv_all)
                                            try:
                                                zall_sum_origin+=float('%.9f'%((fam_npl.ibd_sall-original_dist[2])/original_dist[3]))
                                            except:
                                                pass
                                        else:
                                            zall_sum_origin=z_sum_origin
					for i in range(len(v)):
					    ibd_pair, ibd_all = v[i]
					    if fam_npl.null_mean:
						try:
						    zv_pair = (ibd_pair-fam_npl.null_mean)/fam_npl.null_std
						except ZeroDivisionError:
						    zv_pair = 0
					    else:
						zv_pair = 0
					    zv_pair=float('%.9f'%zv_pair)
					    if pall_flag:
						if fam_npl.sall_null_mean:
						    try:
							zv_all = (ibd_all-fam_npl.sall_null_mean)/fam_npl.sall_null_std
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
                            rtext += 'Z-pairs:{}\tZ-all:{}\n'.format(z_sum_origin/math.sqrt(fam_num),zall_sum_origin/math.sqrt(fam_num))
			    if cz_dist == []:
				if args.exact:              #random sample
                                    if args.verbose>=0:
                                        screen_output.run_out('randomly sampling..')
                                    cz_dist,csall_dist=utils.randomsample(z_dist[idx_apv],args.n_jobs,rep)
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
                            if args.exact:
			        rtext += 'S-pairs:\n'
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
			    if zall_sum_precise:
				Z_all = zall_sum_precise/math.sqrt(fam_num)
				#pall_norm = stats.norm.sf(Z_all)
			    if args.exact:
                                rtext += "S-all:\n"
				rtext += "empiric p-value:"+repr(pr_sall)+"\n"
			    ####################
			pr_s_min = pr_s if pr_s<pr_s_min else pr_s_min
			pr_sall_min = pr_sall if pr_sall<pr_sall_min else pr_sall_min
			possible_min = 1/(int(rep)+1)
			if (pr_s_min,pr_sall_min) == (possible_min,possible_min):
			    break
		    result.write(rtext)
		    result.close()
		    pv = file("{}/pvalue.txt".format(args.output),'a')
		    if args.exact:
			pv.write('%s\t%.9f\t%.9f\n'%(genename,pr_s_min,pr_sall_min))
		    else:
			p_pair_min,p_all_min=min([x for x in asymp_pv if x!=1])
			pv.write('%s\t%.9f\t%.9f\n'%(genename,p_pair_min,p_all_min))
		    pv.close()

