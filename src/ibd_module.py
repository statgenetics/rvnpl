#!/usr/bin/python 
#module containing functions that can calculate the alleles IBD between pairs of 1st cousin, sibpair, half-sibpair and uncle-nephew
#with options to calculate sharing on RV only
from __future__ import division
from RVNPLcpp import ibd_m_cpp
from RVNPLcpp import ibd_rv_cpp

def cousin_ibd(geno,rv_flag=False):    
	#calculate the number of alleles IBD between cousins
	#geno = GT for [grandparents, fam1, fam2] important parent put in the first place of each family
        if rv_flag:
	    ibd = ibd_rv_cpp.cousin_apply(geno)
	else:
	    ibd = ibd_m_cpp.cousin_apply(geno)
        return ibd

def sib_ibd(geno,ibdall=False,rv_flag=False):                  
	#calculate number of alleles IBD between sibpair
	#geno = GT for [father mother sibpair]
	if rv_flag:
	    ibd = ibd_rv_cpp.sib_apply(geno)
	else:
            ibd = ibd_m_cpp.sib_apply(geno)
	if not ibdall:
	    return ibd[1]+ibd[2]*2
	return ibd
def child_ibd(geno):
        #geno = GT for [affected child, affected parent, the other parent]
        pos=[]
        possible_count=0
        rv_shared_count=0
        for allele in geno[:2]:
            tmp_pos=[]
            for idx,parent_a in enumerate(geno[2:]):
                if allele==parent_a:
                    tmp_pos.append(idx)
            pos.append(tmp_pos)
        for p in pos[0]:
            for q in pos[1]:
                if p<2 and q>1:
                    possible_count+=1
                    if geno[0]!=1:
                        rv_shared_count+=1
                if p>1 and q<2:
                    possible_count+=1
                    if geno[1]!=1:
                        rv_shared_count+=1
        return rv_shared_count/possible_count

def hsib_ibd(geno,rv_flag=False):
	#calculate the IBD between half-sibs 
	#geno = GT for [shared parent, 1st other parent, 1st kid, 2nd other parent, 2nd kid]
	ap = geno[:4] 
	ag = geno[4:6]
	bp = geno[:2]+geno[6:8]    
	bg = geno[-2:]
	pg = []
	pg.append(ap)
	pg.append(bp)
	pos_alle = []
	total_p=0
	for g in ag:                                #pick the possible sharing alleles
		if g in bg:
			if rv_flag and g == 1:
			  	continue
			pos_alle.append(g)
	pos_alle = set(pos_alle)
	for a in pos_alle:
		for allele in [m for m in geno[4:6] if m == a]:#
			for alle in [n for n in geno[-2:] if n == a]:   #
				flag = 0
				ag = []+geno[4:6]
				bg = []+geno[-2:]
				f_inherit = []
				for fam in range(2): #for each nuclear family containing affecteds
					possible_inherit=[]
					pos = [i for i,x in enumerate(pg[fam]) if x == a]
					ag.remove(a) if fam==0 else bg.remove(a)
					alt = ag[0] if fam==0 else bg[0]
					pos_alt = [i for i,x in enumerate(pg[fam]) if x == alt]  #the corresponding positions of alleles in parents GT
					for p in pos:
						for q in pos_alt:
							if p<2 and q>1 or p>1 and q<2:
								possible_inherit.append([p,q])    #possible inheritance pattern from parents
					count = 0
					for inherit in possible_inherit:
						if inherit[0] < 2:                  #if the possible sharing allele is from the important parent 
							count +=1
					f_in = count/len(possible_inherit)      #possibility of inheriting sharing allele from important parent
					f_inherit.append(f_in)
					if f_in > 0:              
						flag+=1
				if flag ==2:           # both families have a possibility to inherit the allele
					total_p += f_inherit[0]*f_inherit[1]   #total possibility of half sibs sharing one allele IBD
				else: total_p += 0
	if len(set(geno[:2])) == 1:            # if the share parent is homozygotes, then there is 0.5 chance to share one specific allele
		total_p = total_p/2
	return total_p#/len(pos_alle) if len(pos_alle)>0 else 0

def un_ibd(geno,rv_flag=False):                        
	#estimate IBD for Uncle-Nephew pair
	#geno = GT for [grandparents uncle father mother kid(nephew)]
	if rv_flag:
		ibd=ibd_rv_cpp.un_apply(geno)
	else:
		ibd=ibd_m_cpp.un_apply(geno)
        return ibd

def gp_ibd(geno,rv_flag=False):
        #estimate IBD for grandparent-grandson pair
        #geno = GT for [grandparents father mother kid]
        aff_gp = geno[:2]
        aff_grandson = geno[-2:]
        gp_gt = geno[:4]
        parents_gt = geno[4:8]
        possible_allele=[]
        total_p=0
        for tmp_a in aff_grandson:
            if tmp_a in aff_gp and tmp_a not in possible_allele:
                if rv_flag and tmp_a==1:
                    continue
                possible_allele.append(tmp_a)
        if len(possible_allele)==0:
            return total_p
        for ref in possible_allele:
            #for each possible shared allele
            if ref not in parents_gt:
                continue
            for allele in [tmpa for tmpa in aff_grandson if tmpa==ref]:
                f=0
                m,n=0,0
                bg=[]+aff_grandson
                bg.remove(ref)
                alt=bg[0]
                pos1 = [i for i,x in enumerate(parents_gt) if x == ref]
                pos2 = [i for i,x in enumerate(parents_gt) if x == alt]  #the corresponding positions of alleles in parents GT
                #print pos1,pos2
                for p in pos1:
                        for q in pos2:
                                if p<2 and q>1:
                                        m+=1
                                elif p>1 and q<2:
                                        n+=1
                try:
                    f = m/(m+n)
                except ZeroDivisionError:
                    f = 0
                if f==0:
                    continue
                else:
                    prob=0
                    fm,fn=0,0
                    fbg=[]+parents_gt[:2]
                    fbg.remove(ref)
                    alt=fbg[0]
                    pos1=[i for i,x in enumerate(gp_gt) if x ==ref]
                    pos2=[i for i,x in enumerate(gp_gt) if x ==alt]
                    for p in pos1:
                        for q in pos2:
                            if p<2 and q>1:
                                fm+=1
                            elif p>1 and q<2:
                                fn+=1
                    try:
                        prob=fm/(fm+fn)
                    except ZeroDivisionError:
                        prob=0
                    total_p+=f*prob
        return total_p
