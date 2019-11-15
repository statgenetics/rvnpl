#!/usr/bin/env python2
from __future__ import division
import os
import re
import random
from multiprocessing import Pool
from zipfile import ZipFile
from itertools import repeat
from RVNPL import screen_output
def smart_float(x):
    try:
        return float('%.9f'%x)
    except:
        return None

def smart_int(x):
    try:
        return int(x)
    except ValueError:
        return 0

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

def load_data(fped,fdat,dirpath,args,qtl_flag=False):
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
        #take in ped files
        f = file(pedfile,'r').read()
        lines = f.split('\n')
        lines.pop()
        family = {}
        fam_ids = []
        marker_n=0
        for line in lines:
            tmp=line.split()
            for idx in range(4,len(tmp)):
                if idx==5 and qtl_flag:
                    try:
                        tmp[idx]=float('%.4f'%float(tmp[idx]))
                    except:
                        tmp[idx]=tmp[idx]
                else:
                    tmp[idx]=smart_int(tmp[idx])
            fid = tmp[0]                                   #family id
            if fid not in fam_ids:
                fam_ids.append(fid)
                family[fid]=[]
            marker_n = int((len(tmp)-6)/2)
            #family=[[[fam1member1][fam1member2]...],[fam2member1][]...]
            family[fid].append(tmp[1:])
        #Extract marker freq info from cache
        mafs={}
        archive=ZipFile('{}/cache/{}.cache'.format(pdirname,prefix_name),'r')
        if args.snv:
            try:
                with archive.open('CACHE/{}.{}.freq'.format(prefix_name,chr_id)) as cache_fh:
                    for line in cache_fh:
                        tmp_elements=line.split()
                        mname=tmp_elements[1]
                        tmp_freq=[float(x) for x in tmp_elements[2:]]
                        if mname not in mafs:
                            mafs[mname]=tmp_freq
            except:
                pass
        else:
            #get the freq file for CHP markers
            for mk_name in markername:
                mafs[mk_name]={}
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
                        mafs[tmp[1]][tmp[0]]=CHP_freq
            except:
                pass
            return family, mafs, markername, fam_ids

def ransample_qtl(i):
    numerator=0
    denominator=0
    Q,t=0,0
    if isinstance(null_nd_global[0],dict):
        for fam_dist in null_nd_global:
            r = random.uniform(0,1)
            s = 0
            for item,prob in fam_dist.iteritems():
                s += prob
                if s >= r:
                    n,d=item
                    break
            numerator +=n
            denominator +=d
    else:
        for fam_dist in null_nd_global:
            n,d=fam_dist[random.choice(xrange(len(fam_dist)))]
            numerator +=n
            denominator +=d
    try:
        Q = numerator/denominator
    except:
        Q = 0
    if denominator>0:
        if Q>1:
            t=denominator
        elif Q>0:
            t = Q**2*denominator
    return t

def init_child(data_):
    global null_nd_global
    null_nd_global = data_

def randomsample_qtl(null,n_jobs,rep=20000):
    screen_output.run_out('random sampling...')
    cdist = []
    ##multiple process####
    p = Pool(processes=n_jobs, initializer=init_child, initargs=(null,))
    result = p.imap(ransample_qtl,range(rep), chunksize=rep//args.n_jobs)
    #result.wait()
    p.close()
    p.join()
    cdist = [r for r in result]
    return sorted(cdist)
