#!/usr/bin/env python2 
#combine family class for qualitative and quantitative traits
from __future__ import division
from ctypes import *
import string
import copy
from RVNPLcpp import cmissingparents,cmissing_infer
from RVNPL import screen_output
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
                self.missing_all=self.findall_miss()
                for iid in self.fam_dict.keys():                      #for each fam member
                    if self.fam_dict[iid]['trait']==2 and iid not in self.missing_all:
                            #pick out all affected and sequenced individuals
                            self.affected.append(iid)
                    if self.fam_dict[iid]['parents']!= ['0','0']:
                            self.nonfounder.append(iid)
                    elif self.fam_dict[iid]['parents']== ['0','0']:
                            self.founder.append(iid)
                    self.parents[iid]=self.fam_dict[iid]['parents']
                    self.mates[iid]=self.fam_dict[iid]['mate']
                    self.offspring[iid]=self.fam_dict[iid]['offspring']
                    self.sex[iid]=self.fam_dict[iid]['sex']
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
            self.fallele = []

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

        def dsib(self,aid,bid):                                   #check if it is sibpair
                i = False
                ap=self.fam_dict[aid]['parents']
                bp=self.fam_dict[bid]['parents']
                if ap == bp and ap != ['0','0']:
                        i = True
                return i

        def dhsib(self,aid,bid):
                i = False
                ap=self.fam_dict[aid]['parents']
                bp=self.fam_dict[bid]['parents']
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
                                    elif self.dsib(aid,bid):                              # sibpair
                                            self.sib.append([aid,bid])
                                    elif self.dhsib(aid,bid):
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

class QFamily(Family):
    def __init__(self,mid):
        Family.__init__(self,mid)
        self.hcousin=[]
        self.second_cousin=[]
        self.h2nd_cousin=[]
        self.hun=[]
        self.fcousin_rm=[]
        self.hfcousin_rm=[]
        self.gun=[]
        self.hgp=[]
        self.ggp=[]
        self.offspring_pattern={}
        self.all_members=[]
        self.expect_pair_ibd=[]
        self.pairs=[]
        self.ibdall=[]                  #estimated IBD value based on marker information
        self.prioribd = {}          #list of values of the 1st part in calculate CovTT
        self.prior = []        #sum of prior prob*tau(i,j)*tau(k,l)
        self.post = []        #sum of post prob*tau(i,j)*tau(k,l)
#setdict() is the same
    def iter_offspring(self,nf):
        output=[1 if self.fam_dict[nf]['trait']!='NA' else 0] #indicate trait status
        offspring_with_trait=[x for x in self.fam_dict[nf]['offspring'] if self.fam_dict[x]['trait']!='NA']
        if self.fam_dict[nf]['offspring']:
            mates=self.fam_dict[nf]['mate']
            if len(mates)==1:
                output.append(1 if self.fam_dict[mates[0]]['trait']!='NA' else 0)
                for tmp_off in sorted(self.fam_dict[nf]['offspring'],key=lambda x:len(self.fam_dict[x]['offspring'])):
                    output.append(self.iter_offspring(tmp_off))
                return tuple(output)
            else:
                for mate in mates:
                    tmp_output=[1 if self.fam_dict[mate]['trait']!='NA' else 0]
                    mate_offspring_with_trait=[x for x in self.fam_dict[mate]['offspring'] if self.fam_dict[x]['trait']!='NA']
                    for tmp_off in sorted(self.fam_dict[mate]['offspring'],key=lambda x:len(self.fam_dict[x]['offspring'])):
                        tmp_output.append(self.iter_offspring(tmp_off))
                    output.append(tuple(tmp_output))
                return output
        else:
            if self.fam_dict[nf]['trait']!='NA':
                return (1,)
            else:
                return (0,)

    def core_struct(self,nf):
        #get fam affected structure
        founder_gt=[]
        self.invnf.append(nf)
        self.all_members.append(nf)
        if self.fam_dict[nf]['offspring']:
            mates=self.fam_dict[nf]['mate']
            if len(mates)==1:
                self.all_members.append(mates[0])
                founder_gt.append(sorted(self.fam_dict[mates[0]]['gt'][2*self.mid:2*self.mid+2]))
                for tmp_off in sorted(self.fam_dict[nf]['offspring'],key=lambda x:(len(self.fam_dict[x]['offspring']),self.offspring_pattern[x])):
                    f_gt=self.core_struct(tmp_off)
                    founder_gt.append(f_gt)
            else:
                for mate in mates:
                    mate_fgt=[sorted(self.fam_dict[mate]['gt'][2*self.mid:2*self.mid+2])]
                    self.all_members.append(mate)
                    for tmp_off in sorted(self.fam_dict[mate]['offspring'],key=lambda x:(len(self.fam_dict[x]['offspring']),self.offspring_pattern[x])):
                        f_gt=self.core_struct(tmp_off)
                        mate_fgt.append(f_gt)
                    founder_gt.append(mate_fgt)
        return founder_gt

    def set_famstruct(self):
        # get the affecteds and nonfounder   
        self.famstruct,self.detailstruct,self.parents,self.mates,self.offspring,self.sex,self.offspring_pattern={},{},{},{},{},{},{}
        self.invnf,self.parentsid,self.all_members,self.nonfounder,self.founder,self.missing_all,self.pairs=[],[],[],[],[],[],[]
        for iid in self.fam_dict.keys():                      #for each fam member
            if self.fam_dict[iid]['parents']!= ['0','0']:
                    self.nonfounder.append(iid)
            elif self.fam_dict[iid]['parents']== ['0','0']:
                    self.founder.append(iid)
            self.parents[iid]=self.fam_dict[iid]['parents']
            self.mates[iid]=self.fam_dict[iid]['mate']
            self.offspring[iid]=self.fam_dict[iid]['offspring']
            self.sex[iid]=self.fam_dict[iid]['sex']
        self.missing_all=self.findall_miss()
        #set missing genotype individuals also missing traits
        for iid in self.missing_all:
            self.fam_dict[iid]['trait']='NA'
        #remove missing members that cannot be inferred
        for tmp_p in self.fam_dict.keys():
            if tmp_p in self.missing_all and self.fam_dict[tmp_p]['offspring']==[]:
                self.remove(tmp_p)
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
        if set(allgt)==set([0,1]): #and self.rvibd:
            self.wt_fam = True
        nf_parent={}
        for ind in self.fam_dict.keys():
            nf_p=list(set(self.fam_dict[ind]['parents'])&set(self.nonfounder))
            nf_parent[ind]=nf_p[0] if len(nf_p) else '0'
        self.parentsid = [p for nf in self.nonfounder for p in self.fam_dict[nf]['parents']]
        self.parentsid = sorted(self.parentsid,key=lambda x: self.fam_dict[x]['sex'])
        first_gen=sorted([iid for iid in self.founder if self.generation(iid)==0])
        if min([self.generation(nf) for nf in self.nonfounder]) > -3 and len(first_gen)==2:
            self.simple=True
            self.gpid,self.gmid=sorted(first_gen, key=lambda x: self.fam_dict[x]['sex'])
        self.famstruct['info']=[]
        self.famstruct['offspring']=[]
        for iid in self.fam_dict.keys():
            self.offspring_pattern[iid]=self.iter_offspring(iid)
        if len(first_gen)==2:
            tmp_fgt_info=[]
            self.all_members=first_gen
            self.famstruct['offspring']=self.offspring_pattern[self.all_members[0]]
            for tmp_fg in first_gen:
                tmp_fgt_info.append(sorted(self.fam_dict[tmp_fg]['gt'][2*self.mid:2*self.mid+2]))
            self.famstruct['info']=tmp_fgt_info
            for nf in sorted(self.fam_dict[first_gen[0]]['offspring'],key=lambda x:(len(self.fam_dict[x]['offspring']),self.offspring_pattern[x])):
                founder_info = self.core_struct(nf)
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
            self.famstruct['offspring']=self.offspring_pattern[shared_parent]
            self.all_members.append(shared_parent)
            self.famstruct['info']=[sorted(self.fam_dict[shared_parent]['gt'][2*self.mid:2*self.mid+2])]
            for first_gen_p in first_gen_parents:
                other_parent=first_gen_p[0] if first_gen_p[1]==shared_parent else first_gen_p[1]
                self.all_members.append(other_parent)
                tmp_founder_info=[sorted(self.fam_dict[other_parent]['gt'][2*self.mid:2*self.mid+2])]
                for nf in sorted([iid for iid in self.fam_dict[first_gen_p[0]]['offspring'] if sorted(self.fam_dict[iid]['parents'])==first_gen_p],key=lambda x:(len(self.fam_dict[x]['offspring']),self.offspring_pattern[x])):
                    founder_info = self.core_struct(nf)
                    tmp_founder_info.append(founder_info)
                self.famstruct['info'].append(tmp_founder_info)
        self.invnf_sorted = sorted(self.invnf, key=lambda x: self.generation(x), reverse=True)
        # set value for self.pairs, which should be all pairs of relatives that have phenotype values
        self.all_members=[x for x in self.all_members if self.fam_dict[x]['trait']!='NA']
        all_member_count=len(self.all_members)
        for i in range(all_member_count):
            for j in range(i+1,all_member_count):
                aid,bid=self.all_members[i],self.all_members[j]
                self.pairs.append(sorted([aid,bid]))
        paircount = len(self.pairs)
        if paircount==0:
            self.err=1
        #print self.all_members
        #print self.pairs
        self.expect_pair_ibd=[0 for x in range(paircount)]
        self.prior = [[0 for y in range(paircount)] for x in range(paircount)]
        self.post = [[0 for y in range(paircount)] for x in range(paircount)]

    def relative(self,aid,bid):
        a_gen=self.generation(aid)
        b_gen=self.generation(bid)
        if a_gen==b_gen:
            #possible pairs: sibling, half-sib, cousin, second-cousin
            ap = self.fam_dict[aid]['parents']
            bp = self.fam_dict[bid]['parents']
            if ap!=['0','0'] and bp!=['0','0']:
                if ap==bp:
                    self.sib.append([aid,bid])
                elif set(ap)&set(bp)!=set([]):
                    self.hsib.append([aid,bid])
                else:
                    a_p=[x for x in ap if x in self.nonfounder][0]
                    b_p=[x for x in bp if x in self.nonfounder][0]
                    if a_gen==-2:
                        if self.dsib(a_p,b_p):
                            self.cousin.append([aid,bid])
                        elif self.dhsib(a_p,b_p):
                            self.hcousin.append([aid,bid])
                            self.simple=False
                    elif a_gen==-3:
                        a_gp=[x for x in self.fam_dict[a_p]['parents'] if x in self.nonfounder][0]
                        b_gp=[x for x in self.fam_dict[b_p]['parents'] if x in self.nonfounder][0]
                        if self.dsib(a_gp,b_gp):
                            self.second_cousin.append([aid,bid])
                        elif self.dhsib(a_gp,b_gp):
                            self.h2nd_cousin.append([a_gp,b_gp])
        else:
            senior,junior = (aid,bid) if a_gen > b_gen else (bid,aid)
            senior_gen,junior_gen = (a_gen,b_gen) if a_gen>b_gen else (b_gen,a_gen)
            if senior_gen-junior_gen==1:
                #possible pairs: parent, uncle-nephew, first-cousin once removed or no relation
                jp = self.fam_dict[junior]['parents']
                if jp!=['0','0']:
                    if senior in jp:
                        self.offspring_pairs.append([junior,senior])
                    elif junior_gen!=-1:
                        sp = self.fam_dict[senior]['parents']
                        t_parent=[x for x in jp if x in self.nonfounder][0]
                        if junior_gen==-2 and sp != ['0','0']:
                            if self.dsib(t_parent,senior):
                                self.un.append([senior,junior])
                            elif self.dhsib(t_parent,senior):
                                self.hun.append([senior,junior])
                                self.simple=False
                        elif junior_gen==-3 and sp != ['0','0']:
                            if self.dsib(t_parent,senior):
                                self.un.append([senior,junior])
                            else:
                                gp = [x for x in self.fam_dict[t_parent]['parents'] if x in self.nonfounder][0]
                                tgp = [x for x in sp if x in self.nonfounder][0]
                                if self.dsib(gp,tgp):
                                    self.fcousin_rm.append([senior,junior])
                                elif self.dhsib(gp,tgp):
                                    self.hfcousin_rm.append([senior,junior])
            if senior_gen-junior_gen==2:
                #possible pairs: grandparents, granduncle-nephew, grandhsib-parents
                jp = self.fam_dict[junior]['parents']
                if jp!=['0','0']:
                    t_parent = [x for x in jp if x in self.nonfounder][0]
                    if senior in self.fam_dict[t_parent]['parents']:
                        self.gp.append([senior,junior])
                    else:
                        sp = self.fam_dict[senior]['parents']
                        if sp != ['0','0']:
                            gp = [x for x in self.fam_dict[t_parent]['parents'] if x in self.nonfounder][0]
                            if self.dsib(senior,gp):
                                self.gun.append([senior,junior])
                            elif self.dhsib(senior,gp):
                                self.hgp.append([senior,junior])
            if senior_gen-junior_gen==3:
                self.ggp.append([senior,junior])

        def remove(self,iid):
            Family.remove(self,iid)
            try:
                self.all_members.remove(iid)
            except:
                pass

    def classify_affect(self,m=-1,alleles={},conditional_prob={}):
        #m: marker index
        ibdall = []
        nfnum = len(self.nonfounder)
        #########
        for aid,bid in self.pairs:
            self.relative(aid,bid)
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

    def corre(self,pair,ibdall=[]):
        #return the correlation of traits between relative pairs
        aid, bid = pair
        r = 0
        h = 0.5      #estimated heritability
        if [aid,bid] in self.cousin or [bid,aid] in self.cousin:
            r = h*2*1/16
        elif [aid,bid] in self.sib or [bid, aid] in self.sib:
            r = h*2*1/4
        elif [aid,bid] in self.un or [bid, aid] in self.un:
            r = h*2*1/8
        elif [aid,bid] in self.hsib or [bid, aid] in self.hsib:
            r = h*2*1/8
        elif aid == bid:    #self
            r = 1#h*2*1/2
        elif [aid,bid] in self.gp or [bid, aid] in self.gp:
            r = h*2*1/8
        elif [aid,bid] in self.offspring_pairs or [bid, aid] in self.offspring_pairs:
            #parent-offspring
            r = h*2*1/4
        elif [aid,bid] in self.hcousin or [bid,aid] in self.hcousin:
            r = h*2*1/32
        elif [aid,bid] in self.second_cousin or [bid,aid] in self.second_cousin:
            r = h*2*1/64
        elif [aid,bid] in self.fcousin_rm or [bid,aid] in self.fcousin_rm:
            r = h*2*1/32
        elif [aid,bid] in self.h2nd_cousin or [bid,aid] in self.h2nd_cousin:
            r = h*2*1/128
        elif [aid,bid] in self.hfcousin_rm or [bid,aid] in self.hfcousin_rm:
            r = h*2*1/64
        elif [aid,bid] in self.hun or [bid,aid] in self.hun:
            r = h*2*1/16
        elif [aid,bid] in self.gun or [bid,aid] in self.gun:
            r = h*2*1/16
        elif [aid,bid] in self.hgp or [bid,aid] in self.hgp:
            r = h*2*1/16
        elif [aid,bid] in self.ggp or [bid,aid] in self.ggp:
            r = h*2*1/16
        return r

    def expect(self,pair):       #return the expected IBD between a relative pair
        aid, bid = pair
        tau = 0
        if [aid,bid] in self.cousin or [bid,aid] in self.cousin:
            tau = 0.125
        elif [aid,bid] in self.sib or [bid, aid] in self.sib:
            tau = 0.5
        elif [aid,bid] in self.un or [bid, aid] in self.un:
            tau = 0.25
        elif [aid,bid] in self.hsib or [bid, aid] in self.hsib:
            tau = 0.25
        elif [aid,bid] in self.gp or [bid, aid] in self.gp:
            tau = 0.25
        elif [aid,bid] in self.offspring_pairs or [bid, aid] in self.offspring_pairs: #parent-offspring
            tau = 0.5
        elif [aid,bid] in self.hcousin or [bid,aid] in self.hcousin:
            tau = 1/16
        elif [aid,bid] in self.second_cousin or [bid,aid] in self.second_cousin:
            r = 1/32
        elif [aid,bid] in self.fcousin_rm or [bid,aid] in self.fcousin_rm:
            r = 1/16
        elif [aid,bid] in self.h2nd_cousin or [bid,aid] in self.h2nd_cousin:
            r = 1/64
        elif [aid,bid] in self.hfcousin_rm or [bid,aid] in self.hfcousin_rm:
            r = 1/32
        elif [aid,bid] in self.hun or [bid,aid] in self.hun:
            r = 1/8
        elif [aid,bid] in self.gun or [bid,aid] in self.gun:
            r = 1/8
        elif [aid,bid] in self.hgp or [bid,aid] in self.hgp:
            r = 1/8
        elif [aid,bid] in self.ggp or [bid,aid] in self.ggp:
            r = 1/8
        return tau
