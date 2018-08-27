#include <boost/python.hpp>
#include <boost/python/dict.hpp>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <tuple>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <bitset>
#include <math.h>
#include <assert.h>
//infer missing based on whole family
using namespace std;
typedef std::vector<int> intVec;
typedef std::vector<intVec> VecintVec;
typedef std::vector<VecintVec> VecVecintVec;
typedef std::vector<std::string> strVec;
typedef std::vector<strVec> VecstrVec;
typedef std::vector<long double> doubleVec;
typedef std::pair<strVec,strVec> strVecPair;
typedef map<std::string,string> strMap;
typedef map<std::string,char> charMap;
typedef map<std::string,strVec> strstrMap;
typedef map<std::string,VecintVec> vMap;
typedef map<std::string,long double> strdoubleMap;
typedef map<std::string,strVecPair> strstrVecPairMap;
typedef map<std::string,intVec> alleleMap;
typedef map<std::string,VecVecintVec> myMap;
typedef std::tuple<strdoubleMap,strVec> mytuple;

int convert_to_int(char x)
{
	return (int)x-96;
}

char nth_letter(int n)
{
    assert(n >= 1 && n <= 26);
    return "abcdefghijklmnopqrstuvwxyz"[n-1];
}


strVec conv(strstrMap*a, strVec*b)
{
        strVec result;
        strstrMap tmp_dic=*a;
        strVec sorted_key=*b;
        strVec pred;
        strVec con;
        if(sorted_key.size() > 2){
                string first_key = sorted_key.front();
                strstrMap::iterator it = tmp_dic.find(first_key);
                pred=it->second;
                tmp_dic.erase(it);
                sorted_key.erase(sorted_key.begin());
                con = conv(&tmp_dic,&sorted_key);
        }
        else if(sorted_key.size() == 2){
                strstrMap::iterator it_first = tmp_dic.find(sorted_key.at(0));
                strstrMap::iterator it_second = tmp_dic.find(sorted_key.at(1));
                pred=it_first->second;
                con=it_second->second;
        }
        else if(sorted_key.size() == 1){
                strstrMap::iterator it_first = tmp_dic.find(sorted_key.at(0));
                return it_first->second;
        }
        for (auto iter1=pred.begin(); iter1 != pred.end(); iter1++){
                for (auto iter2=con.begin(); iter2 != con.end(); iter2++){
                        result.push_back(*iter1+*iter2);
                }
        }
        return result;
}

VecstrVec conv2(strstrMap*a){
	strstrMap tmp_dic=*a;
	VecstrVec result;
	strVec sorted_keys;
	for (auto it_k=tmp_dic.begin(); it_k!=tmp_dic.end(); it_k++){
		sorted_keys.push_back(it_k->first);
	}
	sort(sorted_keys.begin(),sorted_keys.end());
	strVec pred;
	VecstrVec con;
	if(sorted_keys.size()==1){
		strVec cond_gt=tmp_dic.at(sorted_keys.at(0));
		for (auto it=cond_gt.begin(); it!=cond_gt.end(); it++){
			result.push_back({*it});
		}
		return result;
	}
	else{
		string first_key=sorted_keys.front();
		strstrMap::iterator it=tmp_dic.find(first_key);
		pred=it->second;
		tmp_dic.erase(it);
		sorted_keys.erase(sorted_keys.begin());
		con = conv2(&tmp_dic);
	}
	for (auto iter1=pred.begin(); iter1!=pred.end(); iter1++){
		for (auto iter2=con.begin(); iter2!=con.end(); iter2++){
			strVec tmp_pred={*iter1};
			strVec tmp_con=*iter2;
			tmp_pred.insert(tmp_pred.end(),tmp_con.begin(),tmp_con.end());
			result.push_back(tmp_pred);
		}
	}
	return result;
}

VecstrVec conv3(strstrMap*a, strVec*arg_sorted_keys){
        strstrMap tmp_dic=*a;
        VecstrVec result;
        strVec sorted_keys=*arg_sorted_keys;
        strVec pred;
        VecstrVec con;
        if(sorted_keys.size()==1){
                strVec cond_gt=tmp_dic.at(sorted_keys.at(0));
                for (auto it=cond_gt.begin(); it!=cond_gt.end(); it++){
                        result.push_back({*it});
                }
                return result;
        }
        else{
                string first_key=sorted_keys.front();
                strstrMap::iterator it=tmp_dic.find(first_key);
                pred=it->second;
                tmp_dic.erase(it);
                sorted_keys.erase(sorted_keys.begin());
                con = conv3(&tmp_dic,&sorted_keys);
        }
        for (auto iter1=pred.begin(); iter1!=pred.end(); iter1++){
                for (auto iter2=con.begin(); iter2!=con.end(); iter2++){
                        strVec tmp_pred={*iter1};
                        strVec tmp_con=*iter2;
                        if (tmp_con.size()>1){
                                string tmp_offgt=(*iter1).substr(0,2);
                                for (auto it_con=tmp_con.begin(); it_con!=tmp_con.end(); it_con++){
                                        tmp_offgt+=(*it_con).substr(0,2);
                                }
                                std::set<char> checker(tmp_offgt.begin(), tmp_offgt.end());
                                if(checker.size()>4){continue;}
                        }
                        tmp_pred.insert(tmp_pred.end(),tmp_con.begin(),tmp_con.end());
                        result.push_back(tmp_pred);
                }
        }
        return result;
}


int postInv(strVec*arginvnf,strVec*argfounder,strstrMap*argparents,alleleMap*argalleles){
        strVec invnf=*arginvnf;
        strVec founder=*argfounder;
        strstrMap parents=*argparents;
        alleleMap alleles=*argalleles;
	myMap tmpv;
	vMap v;
	vMap pos;
	int offspring=0;
	strVec branch_order;
	int count=1;
	for (auto it_nf=invnf.begin(); it_nf!=invnf.end(); it_nf++){
		//for each nonfounder
		strVec tmp_parents=parents.at(*it_nf);
		//cout<<"nf:"<<*it_nf<<" "<<endl;
		intVec pgeno;   //parents genotype
		for (auto it_parent=tmp_parents.begin(); it_parent!=tmp_parents.end(); it_parent++){
			intVec t_a=alleles.at(*it_parent);
			pgeno.insert(pgeno.end(),t_a.begin(),t_a.end());
		}
		VecintVec tmpos;
		intVec nf_a=alleles.at(*it_nf);
		for (auto it_nf_a=nf_a.begin(); it_nf_a!=nf_a.end(); it_nf_a++){
			intVec a_pos;
			for (auto it_parent_a=pgeno.begin(); it_parent_a!=pgeno.end(); it_parent_a++){
				if(*it_parent_a==*it_nf_a){a_pos.push_back(it_parent_a-pgeno.begin());}
			}
			tmpos.push_back(a_pos);
		}
		for (auto it_pos1=tmpos.at(0).begin(); it_pos1!=tmpos.at(0).end(); it_pos1++){
			for (auto it_pos2=tmpos.at(1).begin(); it_pos2!=tmpos.at(1).end(); it_pos2++){
				if (*it_pos1<2&&*it_pos2>1||*it_pos1>1&&*it_pos2<2){
					intVec t_pos={*it_pos1,*it_pos2};
					intVec t_pos2={*it_pos2,*it_pos1};
					VecintVec nf_pos;
					if (pos.find(*it_nf) != pos.end()){
						nf_pos = pos.at(*it_nf);
						if (find(nf_pos.begin(),nf_pos.end(),t_pos) == nf_pos.end() &&
						    find(nf_pos.begin(),nf_pos.end(),t_pos2) == nf_pos.end()){
							pos.at(*it_nf).push_back(t_pos);
							}
					}
					else{
						nf_pos.push_back(t_pos);
						pos.insert(std::pair<string,VecintVec>(*it_nf,nf_pos));}
				}
			}
		}
		if (pos.find(*it_nf)==pos.end()){
			return 0;
		}
		count*=pos.at(*it_nf).size();
	}
	return count;
}

strVec nuclear_infer(doubleVec *argmaf,strVec *arg_offspring_GT,strVec *arg_parent_GT){
        strVec possible_GT;
        doubleVec maf = *argmaf;
        strVec offspring_GT = *arg_offspring_GT;
        strVec parent_GT = *arg_parent_GT;
        if (std::find(offspring_GT.begin(),offspring_GT.end(),"00") != offspring_GT.end()){
                //if "00" is found in offspring_GT, then the fam is uninformative
                possible_GT={"0000"};
                return possible_GT;
        }
        else{
                strVec parental_GT;
                strVec alleles;
                for (int i=0; i!=maf.size(); i++){
                        char letter=nth_letter(i+1);
                        string s(1,letter);
                        alleles.push_back(s);
                }
                for (int ait1=0; ait1!=maf.size(); ait1++){
                        for (int ait2=0; ait2!=maf.size(); ait2++){
                                string gt=alleles.at(ait1)+alleles.at(ait2);
                                std::sort(gt.begin(),gt.end());
                                if (find(parental_GT.begin(),parental_GT.end(),gt)==parental_GT.end()){
                                        parental_GT.push_back(gt);
                                }
                        }
                }
                long double prior_offgt = 0;
                strMap combined_prob;
                string fgt=parent_GT.at(0);
                string mgt=parent_GT.at(1);
                strVec possible_fgt, possible_mgt;
                if (parent_GT.at(0) == "None"){
                        possible_fgt = parental_GT;
                }
                else{possible_fgt.push_back(fgt);}
                if (parent_GT.at(1) == "None"){
                        possible_mgt = parental_GT;
                }
                else{possible_mgt.push_back(mgt);}
                for (auto fit=possible_fgt.begin();fit!=possible_fgt.end();fit++){
                        for (auto mit=possible_mgt.begin(); mit!=possible_mgt.end(); mit++){
                                string combined_gt = *fit+*mit;
                                long double conditional_p=1;
                                for (auto off_it=offspring_GT.begin();off_it!=offspring_GT.end();off_it++){
                                        double success=0;
                                        std::vector<int> pos1,pos2;
                                        string tmp_offgt = *off_it;
                                        for (auto tmp_it=combined_gt.begin();tmp_it!=combined_gt.end();tmp_it++){
                                                int idx=tmp_it-combined_gt.begin();
                                                if (*tmp_it==tmp_offgt.at(0)){
                                                        pos1.push_back(idx);
                                                }
                                                if (*tmp_it==tmp_offgt.at(1)){
                                                        pos2.push_back(idx);
                                                }
                                        }
                                        for (auto pit = pos1.begin();pit!=pos1.end();pit++){
                                                for (auto qit = pos2.begin();qit!=pos2.end();qit++){
                                                        if (*pit<2&&*qit>1 || *pit>1&&*qit<2){
                                                                success+=1;
                                                        }
                                                }
                                        }
                                        conditional_p*=success/8;
                                        if (conditional_p==0){break;}
                                }
                                if (conditional_p==0){continue;}
                                possible_GT.push_back(combined_gt);
                        }
                }
                return possible_GT;
        }
}



strVecPair hsib_infer(strstrMap *argcond, string *argshared, int *argshared_pos, strstrMap *arginfer_persons){
	strstrMap cond_dic=*argcond;
	string shared=*argshared;
	strstrMap infer_persons=*arginfer_persons;
	int shared_pos=*argshared_pos;
	strVecPair output;
	strVec mates;
	for (auto it_mate=cond_dic.begin(); it_mate!=cond_dic.end(); it_mate++){
		mates.push_back(it_mate->first);
	}
	sort(mates.begin(),mates.end());
	VecstrVec pre_combined_gt=conv2(&cond_dic);	
	VecstrVec new_combined_gt;
	//filter out those that have conflict genotypes on shared partner
	for (auto it_gt=pre_combined_gt.begin(); it_gt!=pre_combined_gt.end(); it_gt++){
		strVec shared_gt;
		bool conflict_flag=false;
		for (auto it_gt_pair=(*it_gt).begin(); it_gt_pair!=(*it_gt).end(); it_gt_pair++){
			string pair_gt;
			//cout<<*it_gt_pair<<" ";
			int t_len=(*it_gt_pair).length();
			if (t_len>4){
				pair_gt=(*it_gt_pair).substr(t_len-4,4);
			}
			else{
				pair_gt=*it_gt_pair;
			}
			string tmp_shared_gt=pair_gt.substr(2*shared_pos,2);
			//cout<<"shared_gt:"<<tmp_shared_gt<<" ";
			if (find(shared_gt.begin(),shared_gt.end(),tmp_shared_gt)==shared_gt.end() && shared_gt.size()>0){conflict_flag=true;break;}
			else{shared_gt.push_back(tmp_shared_gt);}
		}
		//cout<<endl;
		if (!conflict_flag){
			new_combined_gt.push_back(*it_gt);
		}
	}
	strVec possible_GT;
	for (auto it_gt=new_combined_gt.begin(); it_gt!=new_combined_gt.end(); it_gt++){
		int tmp_len=(*it_gt).at(0).length();
		string tmp_pgt=(*it_gt).at(0).substr(tmp_len-4,4);
		string key=tmp_pgt.substr(2*shared_pos,2);
		for (auto it_tmp_gt=(*it_gt).begin(); it_tmp_gt!=(*it_gt).end(); it_tmp_gt++){
			int t_len = (*it_tmp_gt).length();
			string parent_gt=(*it_tmp_gt).substr(t_len-4,4);
			string mate_gt=parent_gt.substr(2*(1-shared_pos),2);
			key+=mate_gt;
			if (t_len>4){
				key+=(*it_tmp_gt).substr(0,t_len-4);
			}
		}
		possible_GT.push_back(key);
	}
	strVec combined_id={shared};
	for (auto it_mate=mates.begin(); it_mate!=mates.end(); it_mate++){
		combined_id.push_back(*it_mate);
		if (infer_persons.at(*it_mate).size()>2){
			for (auto it=infer_persons.at(*it_mate).begin(); it!=infer_persons.at(*it_mate).end()-2; it++){
				combined_id.push_back(*it);
			}
		}
	}
	output=std::make_pair(possible_GT,combined_id);
	return output;
}

bool check_sufficient(char allele,string *argnonfounder_gt, strVec *argmp_offgt){
        //check if the allele is sufficient
        string nonfounder_gt=*argnonfounder_gt;
        strVec mp_offgt=*argmp_offgt;
        string tmp_gt1, tmp_gt2;
        tmp_gt1+=allele;
        tmp_gt1+=nonfounder_gt.at(0);
        tmp_gt2+=allele;
        tmp_gt2+=nonfounder_gt.at(1);
        std::sort(tmp_gt1.begin(),tmp_gt1.end());
        std::sort(tmp_gt2.begin(),tmp_gt2.end());
        for (auto it_offgt=mp_offgt.begin(); it_offgt!=mp_offgt.end(); it_offgt++){
                string offgt=*it_offgt;
                std::sort(offgt.begin(),offgt.end());
                if (offgt!=tmp_gt1 && offgt!=tmp_gt2){return false;}
        }
        return true;
}

strVecPair infer_engine(doubleVec *argmaf, strVec *argcouple,intVec *argcouple_founder, strdoubleMap *argpossible_null_GT, strVec *argoffspring, strMap *argtmp_alleles, strstrVecPairMap *argpossible_missing_GT, int *argsimple){
	strVec couple=*argcouple;
        intVec couple_founder=*argcouple_founder;
        strdoubleMap possible_null_GT=*argpossible_null_GT;
	strVec offspring=*argoffspring;
	std::sort(offspring.begin(),offspring.end());
	strMap tmp_alleles=*argtmp_alleles;
	doubleVec maf=*argmaf;
	strstrVecPairMap possible_missing_GT=*argpossible_missing_GT;
        int simple=*argsimple;
	strVecPair output;
	strstrMap nuclear_possgt;
	string father=couple.at(0);
	string mother=couple.at(1);	
	strVec mp_pgt;
        bool both_founder=(couple_founder.at(0)==1&&couple_founder.at(1)==1)?true:false;
        int founder_idx;
        if (!both_founder){founder_idx=(couple_founder.at(0)==1)?0:1;}
	if (tmp_alleles.at(father)=="NULL"){
		mp_pgt.push_back("None");
	}
	else{
		mp_pgt.push_back(tmp_alleles.at(father));
	}
	if (tmp_alleles.at(mother)=="NULL"){
		mp_pgt.push_back("None");
	}
	else{
		mp_pgt.push_back(tmp_alleles.at(mother));
	}
	/*cout<<"mp_pgt: ";
	for (auto it=mp_pgt.begin(); it!=mp_pgt.end(); it++){
		cout<<*it;
	}
	cout<<endl;*/
	strstrMap miss_off_gt;
	strVec combined_id;
	for (auto it_off=offspring.begin(); it_off!=offspring.end(); it_off++){
		if (possible_missing_GT.find(*it_off)!=possible_missing_GT.end()){
                        strVec tmp_gt=possible_missing_GT.at(*it_off).first;
                        strVec tmp_ids=possible_missing_GT.at(*it_off).second;
                        int off_id=find(tmp_ids.begin(),tmp_ids.end(),*it_off)-tmp_ids.begin();
                        //put self at first
                        strVec new_gt;
                        strVec new_ids={*it_off};
                        for (int it_id=0; it_id!=tmp_ids.size(); it_id++){
                                if(it_id!=off_id){new_ids.push_back(tmp_ids.at(it_id));}
                        }
                        for (auto itt=tmp_gt.begin(); itt!=tmp_gt.end(); itt++){
                                string tmp_new_gt=(*itt).substr(2*off_id,2);
                                for (int itt_id=0; itt_id!=tmp_ids.size(); itt_id++){
                                        if (itt_id!=off_id){tmp_new_gt+=(*itt).substr(2*itt_id,2);}
                                }
                                new_gt.push_back(tmp_new_gt);
                        }
                        miss_off_gt.insert(std::pair<string,strVec>(*it_off,new_gt));
                        combined_id.insert(combined_id.end(),new_ids.begin(),new_ids.end());
		}
		else{
			miss_off_gt.insert(std::pair<string,strVec>(*it_off,{tmp_alleles.at(*it_off)}));
			combined_id.push_back(*it_off);
		}
	}
	/*cout<<"miss_off_gt: ";
	for (auto it_k=miss_off_gt.begin(); it_k!=miss_off_gt.end(); it_k++){
		for (auto itt=it_k->second.begin(); itt!=it_k->second.end(); itt++){
			cout<<*itt;
		}
		cout<<" ";
	}
	cout<<endl;*/
	VecstrVec return_gt=conv3(&miss_off_gt,&offspring);
	strVec combined_possible_GT;
        intVec combined_id_offspring;
        for (auto it_off=offspring.begin(); it_off!=offspring.end(); it_off++){
                int combined_idx=find(combined_id.begin(),combined_id.end(),*it_off)-combined_id.begin();
                combined_id_offspring.push_back(combined_idx);
        }
	for (auto it_rgt=return_gt.begin(); it_rgt!=return_gt.end(); it_rgt++){
		strVec mp_offgt;
		string mp_offgt_key;
		strVec tmp_possible_GT;
		string str_rgt;
                for (auto itt_rgt=(*it_rgt).begin(); itt_rgt!=(*it_rgt).end(); itt_rgt++){
                        str_rgt+=*itt_rgt;
                }
                for (int it_off_idx=0; it_off_idx!=offspring.size(); it_off_idx++){
                        string tmp_offgt=str_rgt.substr(2*combined_id_offspring.at(it_off_idx),2);
			mp_offgt.push_back(tmp_offgt);
			mp_offgt_key+=tmp_offgt;
		}
                if (nuclear_possgt.find(mp_offgt_key)!=nuclear_possgt.end()){
                        tmp_possible_GT=nuclear_possgt.at(mp_offgt_key);
                }
                else{
                        tmp_possible_GT=nuclear_infer(&maf,&mp_offgt,&mp_pgt);
                        nuclear_possgt.insert(std::pair<string,strVec>(mp_offgt_key,tmp_possible_GT));
                }
                if (!tmp_possible_GT.size()){continue;}
                for (auto it_pgt=tmp_possible_GT.begin(); it_pgt!=tmp_possible_GT.end(); it_pgt++){
                        if (!both_founder && simple){
                                //NEED REVISION, ADD HERE: remove founder genotype that have too many RVs
                                string founder_gt;
                                string nonfounder_gt;
                                founder_gt=(*it_pgt).substr(2*founder_idx,2);
                                nonfounder_gt=(*it_pgt).substr(2*(1-founder_idx),2);
				if (simple==1){
					if (check_sufficient(founder_gt.at(0),&nonfounder_gt,&mp_offgt) && maf.at(convert_to_int(founder_gt.at(1))-1)<1E-4 ||
						check_sufficient(founder_gt.at(1),&nonfounder_gt,&mp_offgt) && maf.at(convert_to_int(founder_gt.at(0))-1)<1E-4){
						continue;
					}
				}
				else if (simple==2){
					if (check_sufficient(founder_gt.at(0),&nonfounder_gt,&mp_offgt) && founder_gt.at(1)!='a' ||
						check_sufficient(founder_gt.at(1),&nonfounder_gt,&mp_offgt) && founder_gt.at(0)!='a'){
						continue;
					}
				}
                        }
                        string combined_gt=str_rgt+*it_pgt;
                        combined_possible_GT.push_back(combined_gt);
                }

	}
	combined_id.push_back(father);
	combined_id.push_back(mother);
	output=std::make_pair(combined_possible_GT,combined_id);
	return output;
}

strVecPair initial_infer(doubleVec *argmaf, strVec *argfounder, strdoubleMap *argpossible_null_GT, strstrMap *parents, strstrMap *argmates, strstrMap *argoffspring, charMap *arg_sex, strVec *sorted_miss_persons, strMap *t_alleles, int *argsimple){
	//initial infer on missing GT based on known GT
	strstrMap tparents=*parents;
	strstrMap tmates=*argmates;
	strstrMap toffspring=*argoffspring;
	strVec tsorted_miss_persons=*sorted_miss_persons;
	strMap tmp_alleles=*t_alleles;
	charMap sex=*arg_sex;
	doubleVec maf=*argmaf;
        strVec founder=*argfounder;
        strdoubleMap possible_null_GT=*argpossible_null_GT;
        int simple=*argsimple;
	strstrVecPairMap possible_missing_GT;
	strVecPair output;
	for (auto it_miss=tsorted_miss_persons.begin(); it_miss!=tsorted_miss_persons.end(); it_miss++){
		strVecPair tmp_pair;// ({"NN"},{"NN"});
		if (possible_missing_GT.find(*it_miss)!=possible_missing_GT.end()){continue;}
		strVec mates=tmates.at(*it_miss);
		if (mates.size()==1 && tmates.at(mates.at(0)).size()==1){
			//single marriage
			string mate=mates.at(0);
			if (possible_missing_GT.find(mate)!=possible_missing_GT.end()){
				tmp_pair=possible_missing_GT.at(mate);
				possible_missing_GT.insert(std::pair<string,strVecPair>(*it_miss,tmp_pair));
			}
			else{
				strVec couple={*it_miss,mate};
                                intVec couple_founder;
                                for (auto it_couple=couple.begin(); it_couple!=couple.end(); it_couple++){
                                        int in_founder=find(founder.begin(),founder.end(),*it_couple)==founder.end()?0:1;
                                        couple_founder.push_back(in_founder);
                                }
				strVec offspring=toffspring.at(*it_miss);
				tmp_pair=infer_engine(&maf,&couple,&couple_founder, &possible_null_GT,&offspring,&tmp_alleles,&possible_missing_GT,&simple);
				possible_missing_GT.insert(std::pair<string,strVecPair>(*it_miss,tmp_pair));
				possible_missing_GT.insert(std::pair<string,strVecPair>(mate,tmp_pair));
			}
		}
		else if (mates.size()>1){
			//more than 1 marriage
			strstrMap tmp_cond;
			strstrMap infer_persons;
			for (auto it_mate=mates.begin(); it_mate!=mates.end(); it_mate++){
				strVec couple={*it_miss,*it_mate};
                                intVec couple_founder;
                                for (auto it_couple=couple.begin(); it_couple!=couple.end(); it_couple++){
                                        int in_founder=find(founder.begin(),founder.end(),*it_couple)==founder.end()?0:1;
                                        couple_founder.push_back(in_founder);
                                }
				strVec offspring=toffspring.at(*it_mate);
				tmp_pair=infer_engine(&maf,&couple,&couple_founder, &possible_null_GT,&offspring,&tmp_alleles,&possible_missing_GT,&simple);
				tmp_cond.insert(std::pair<string,strVec>(*it_mate,tmp_pair.first));
				infer_persons.insert(std::pair<string,strVec>(*it_mate,tmp_pair.second));
			}
			int shared_pos = (sex.at(*it_miss)=='1') ? 0 : 1;
			tmp_pair=hsib_infer(&tmp_cond,&(*it_miss),&shared_pos,&infer_persons);
			possible_missing_GT.insert(std::pair<string,strVecPair>(*it_miss,tmp_pair));
			for (auto it_mate=mates.begin(); it_mate!=mates.end(); it_mate++){
				possible_missing_GT.insert(std::pair<string,strVecPair>(*it_mate,tmp_pair));
			}
		}	
	}
	string last_miss_person=tsorted_miss_persons.back();
	possible_missing_GT.insert(std::pair<string,strVecPair>("~combined",possible_missing_GT.at(last_miss_person)));
	for (auto it_miss=tsorted_miss_persons.begin(); it_miss!=tsorted_miss_persons.end(); it_miss++){
		strVec poss_persons=possible_missing_GT.at("~combined").second;
		if (find(poss_persons.begin(),poss_persons.end(),*it_miss)==poss_persons.end() && possible_missing_GT.find(*it_miss)!=possible_missing_GT.end())
		{
			//if some missing person not included in the least confirmed person's GT configurations
			strstrMap combined_prob;
			strVec gt_id;
			strVec possible_gt;
			combined_prob.insert(std::pair<string,strVec>(*it_miss,possible_missing_GT.at(*it_miss).first));
			combined_prob.insert(std::pair<string,strVec>("~combined",possible_missing_GT.at("~combined").first));
			VecstrVec return_gt;
			return_gt=conv2(&combined_prob);
			gt_id.insert(gt_id.end(),possible_missing_GT.at(*it_miss).second.begin(),possible_missing_GT.at(*it_miss).second.end());
			gt_id.insert(gt_id.end(),possible_missing_GT.at("~combined").second.begin(),possible_missing_GT.at("~combined").second.end());
			for (auto it_gt=return_gt.begin(); it_gt!=return_gt.end(); it_gt++){
				string tmp_gt=(*it_gt).at(0)+(*it_gt).at(1);
				possible_gt.push_back(tmp_gt);
			}
			strVecPair tmp_gt_pair (possible_gt,gt_id);
			possible_missing_GT.at("~combined")=tmp_gt_pair;
		}
	}
	
	intVec miss_ids;
	strVec miss_possible_gt;
	strVec miss_gt_ids;
	strVec tmp_gt_id=possible_missing_GT.at("~combined").second;
	strVec tmp_poss_gt=possible_missing_GT.at("~combined").first;
	for (auto it_gt_id=tmp_gt_id.begin(); it_gt_id!=tmp_gt_id.end(); it_gt_id++){
		if (find(tsorted_miss_persons.begin(), tsorted_miss_persons.end(), *it_gt_id)!=tsorted_miss_persons.end()){
			miss_gt_ids.push_back(*it_gt_id);
			miss_ids.push_back(it_gt_id-tmp_gt_id.begin());
		}
	}

	for (auto it_poss_gt=tmp_poss_gt.begin(); it_poss_gt!=tmp_poss_gt.end(); it_poss_gt++){
		string tmp_miss_gt;
		for (auto it_mid=miss_ids.begin(); it_mid!=miss_ids.end(); it_mid++){
			tmp_miss_gt+=(*it_poss_gt).substr(2*(*it_mid),2);
		}
		miss_possible_gt.push_back(tmp_miss_gt);
	}
	output=std::make_pair(miss_possible_gt,miss_gt_ids);
	//cout<<"output size:"<<miss_possible_gt.size()<<endl;
	return output;
}


mytuple probability(doubleVec *argmaf,strVec*arginvnf,strVec*argfounder,strstrMap*argparents,strstrMap*argmates,strstrMap*argoffspring,charMap*argsex,strVec*arg_sorted_miss_persons,alleleMap*argalleles,int *argsimple){
        strVec invnf=*arginvnf;
        strVec founder=*argfounder;
        strstrMap parents=*argparents;
        strstrMap mates=*argmates;
        strstrMap offspring=*argoffspring;
	charMap sex=*argsex;
	strVec sorted_miss_persons=*arg_sorted_miss_persons;
        alleleMap talleles=*argalleles;
        int simple=*argsimple;
        strdoubleMap conditional_prob;
        strdoubleMap tmp_conditional_prob;
	mytuple output;
        doubleVec maf = *argmaf;
	strdoubleMap possible_null_GT;
	strVec alleles;
	strVec null_gt;
	strMap founder_alleles;
	strMap t_alleles;
	string known_gt;
	for (auto it_ind=talleles.begin(); it_ind!=talleles.end(); it_ind++){
		if (std::find(it_ind->second.begin(),it_ind->second.end(),0)==it_ind->second.end()){
			//known genotypes
			char letter=nth_letter(it_ind->second.at(0));
			string tmps(1,letter);
			char letter1=nth_letter(it_ind->second.at(1));
			string tmps1(1,letter1);
			string gt=tmps+tmps1;
			std::sort(gt.begin(),gt.end());
			known_gt+=gt;
			t_alleles.insert(std::pair<string,string>(it_ind->first,gt));
			strVec::iterator it_founder=find(founder.begin(),founder.end(),it_ind->first);
			if(it_founder!=founder.end()){
				founder_alleles.insert(std::pair<string,string>(it_ind->first,gt));
			}
		}
		else{
			strVec::iterator it_founder=find(founder.begin(),founder.end(),it_ind->first);
			t_alleles.insert(std::pair<string,string>(it_ind->first,"NULL"));
			if(it_founder!=founder.end()){
				founder_alleles.insert(std::pair<string,string>(it_ind->first,"NULL"));
			}
		}
	}
	//all possible genotypes for one individual
	for (int i=0; i!=maf.size(); i++){
		char letter=nth_letter(i+1);
		string s(1,letter);
		alleles.push_back(s);
	}
	for (int ait1=0; ait1!=maf.size(); ait1++){
		for (int ait2=0; ait2!=maf.size(); ait2++){
			string gt=alleles.at(ait1)+alleles.at(ait2);
			std::sort(gt.begin(),gt.end());
			if (find(null_gt.begin(),null_gt.end(),gt)!=null_gt.end()){
				long double tmp_maf=maf.at(ait1)*maf.at(ait2);
				possible_null_GT.at(gt)+=tmp_maf;
			}
			else{
				null_gt.push_back(gt);
				long double tmp_maf=maf.at(ait1)*maf.at(ait2);
				possible_null_GT.insert(std::pair<string,long double>(gt,tmp_maf));
			}
		}
	}
	//initial infer on possible GT
	strVecPair possible_GT_pair=initial_infer(&maf, &founder, &possible_null_GT, &parents, &mates, &offspring, &sex, &sorted_miss_persons, &t_alleles, &simple);
	strVec possible_GT=possible_GT_pair.first;
	strVec missing_person_id=possible_GT_pair.second;
	int null_inv_size=pow(2,2*invnf.size());
	for (auto it_poss_gt=possible_GT.begin(); it_poss_gt!=possible_GT.end(); it_poss_gt++){
		// for each possible missing GT
		long double founder_gt_prob=1.0;
		string miss_gt=*it_poss_gt;
		string tmp_all_gt=miss_gt+known_gt;
		//cout<<"For each:"<<endl;
		//cout<<miss_gt<<endl;
		for (int miss_idx=0; miss_idx!=missing_person_id.size(); miss_idx++){
			string tmp_gt=miss_gt.substr(2*miss_idx,2);
			string person_id=missing_person_id.at(miss_idx);
			intVec int_gt={convert_to_int(tmp_gt.at(0)),convert_to_int(tmp_gt.at(1))};
			talleles.at(person_id)=int_gt;
			strVec::iterator it_founder=find(founder.begin(),founder.end(),person_id);
			if(it_founder!=founder.end()){
				founder_alleles.at(person_id)=tmp_gt;
				founder_gt_prob*=possible_null_GT.at(tmp_gt);
			}
		}
		int postinv_size;
 		if (tmp_all_gt.find_first_not_of(tmp_all_gt[0]) == std::string::npos){
 			//uninformative
 			postinv_size=null_inv_size;
 		}
 		else{
 			postinv_size=postInv(&invnf,&founder,&parents,&talleles);
 		}
		//cout<<"  prob: "<<founder_gt_prob<<endl;
		//cout<<"postinv:"<<postinv.size()<<"null_inv:"<<null_inv_size<<endl;
		long double post_prob = (long double) postinv_size/null_inv_size;
		long double combined_prob=founder_gt_prob*post_prob;
		//cout<<"combined_prob:"<<combined_prob<<endl;
		if (combined_prob>0){	
			tmp_conditional_prob.insert(std::pair<string,long double>(*it_poss_gt,combined_prob));
		}
	}
	long double sum=0;
	for (auto it_gt=tmp_conditional_prob.begin(); it_gt!=tmp_conditional_prob.end(); it_gt++){
		sum+=it_gt->second;
	}
	//cout<<"sum:"<<sum<<" count:"<<tmp_conditional_prob.size()<<endl;
	for (auto it_gt=tmp_conditional_prob.begin(); it_gt!=tmp_conditional_prob.end(); it_gt++){
		conditional_prob.insert(std::pair<string,long double>(it_gt->first,it_gt->second/sum));
		//cout<<it_gt->first<<":"<<conditional_prob.at(it_gt->first)<<endl;
	}
	output=std::make_tuple(conditional_prob,missing_person_id);
	return output;
}

/*int main(){
	strVecPair output;
	strstrMap input={{"first",{"abcd","acca","adac"}},{"second",{"abbe","acbf","adbg"}}};
	strstrMap input2={{"first",{"share","first"}},{"second",{"shared","second"}}};
	string shared="share";
	int shared_pos=0;
	strVec couple={"1","2"};
	strMap alleles={{"1","NULL"},{"2","NULL"},{"3","ab"},{"4","ab"}};
	strstrVecPairMap possible_missing_GT;
	doubleVec maf={0.99, 0.01};
	strstrMap parents={{"1",{"0","0"}},{"2",{"0","0"}},{"3",{"1","2"}},{"4",{"1","2"}}};
	strVec sorted_miss_persons={"1","2"};
	strstrMap mates={{"1",{"2"}},{"2",{"1"}}};
	strstrMap offspring={{"1",{"3","4"}},{"2",{"3","4"}}};
	strMap sex={{"1",{"1"}},{"2",{"2"}},{"3",{"1"}},{"4",{"2"}}};
	output=initial_infer(&maf, &parents, &mates, &offspring, &sex, &sorted_miss_persons, &alleles);
	mytuple prob=probability(&maf,&invnf,&founder,&parents,&mates,&offspring,&sex, &sorted_miss_persons, &alleles);
	for (auto it=output.first.begin(); it!=output.first.end(); it++){
		cout<<*it<<endl;
	}
	for (auto it=output.second.begin(); it!=output.second.end(); it++){
		cout<<*it<<" ";
	}
	cout<<endl;
        doubleVec maf={0.9945, 0.0055};
	//doubleVec maf={0.995411905826, 0.000697276427577, 0.00119592940227, 0.00269488834426};
	strVec invnf={"1798", "1796", "1805", "1801", "1813", "1802", "1814", "1804", "1818", "1819", "1820"};
	strVec founder={"1830", "1829", "1797", "1808", "1807", "1803"};
	strstrMap parents={{"1814",{"1803", "1802"}}, {"1807",{"0", "0"}}, {"1805",{"1829", "1830"}}, {"1798",{"1829", "1830"}}, {"1829",{"0", "0"}}, {"1802",{"1829", "1830"}}, {"1801",{"1829", "1830"}}, {"1813",{"1808", "1801"}}, {"1797",{"0", "0"}}, {"1818",{"1807", "1804"}}, {"1830",{"0", "0"}}, {"1796",{"1797", "1798"}}, {"1820",{"1807", "1804"}}, {"1804",{"1829", "1830"}}, {"1803",{"0", "0"}}, {"1808",{"0", "0"}}, {"1819",{"1807", "1804"}}};
	strstrMap mates={{"1807",{"1804"}}, {"1798",{"1797"}}, {"1829",{"1830"}}, {"1802",{"1803"}}, {"1801",{"1808"}}, {"1797",{"1798"}}, {"1830",{"1829"}}, {"1804",{"1807"}}, {"1803",{"1802"}}, {"1808",{"1801"}}};
	strstrMap offspring={{"1804",{"1818","1819","1820"}}, {"1807",{"1818","1819","1820"}}, {"1798",{"1796"}}, {"1797",{"1796"}}, {"1829",{"1798","1805","1801","1802","1804"}}, {"1830",{"1798","1805","1801","1802","1804"}}, {"1802",{"1814"}}, {"1803",{"1814"}}, {"1801",{"1813"}}, {"1808",{"1813"}}};
	strMap sex={{"1814",{"1"}}, {"1807",{"1"}}, {"1805",{"2"}}, {"1798",{"2"}}, {"1829",{"1"}}, {"1802",{"2"}}, {"1801",{"2"}}, {"1813",{"2"}}, {"1797",{"1"}}, {"1818",{"2"}}, {"1830",{"2"}}, {"1796",{"1"}}, {"1820",{"2"}}, {"1804",{"2"}}, {"1803",{"1"}}, {"1808",{"1"}}, {"1819",{"2"}}};
	strVec sorted_miss_persons={"1797","1801","1808","1802","1803","1804","1807","1829","1830"};
	alleleMap alleles={{"1814",{1, 2}}, {"1807",{0, 0}}, {"1805",{1, 1}}, {"1798",{1, 1}}, {"1829",{0, 0}}, {"1802",{0, 0}}, {"1801",{0, 0}}, {"1813",{1, 1}}, {"1797",{0, 0}}, {"1818",{1, 1}}, {"1830",{0, 0}}, {"1796",{1, 1}}, {"1820",{1, 1}}, {"1804",{0, 0}}, {"1803",{0, 0}}, {"1808",{0, 0}}, {"1819",{1, 1}}};
	strdoubleMap prob=probability(&maf,&invnf,&founder,&parents,&mates,&offspring,&sex, &sorted_miss_persons, &alleles);
}*/

doubleVec todoubleVec(boost::python::list pylist){
        doubleVec clist;
        for (int i=0; i<len(pylist);i++){
                long double gt=boost::python::extract<double>(pylist[i]);
                clist.push_back(gt);
        }
        return clist;
}

strVec tostrVec(boost::python::list mylist){
        strVec output;
        for (int i=0; i<len(mylist); i++){
                string value=boost::python::extract<string>(mylist[i]);
                output.push_back(value);
        }
        return output;
}

intVec tointVec(boost::python::list mylist){
        intVec output;
        for (int i=0; i<len(mylist); i++){
                int value=boost::python::extract<int>(mylist[i]);
                output.push_back(value);
        }
        return output;
}

strstrMap tomap(boost::python::dict mydict){
        strstrMap convert_map;
        boost::python::list keys = mydict.keys();
        for(int i = 0; i<len(keys);++i){
                std::string key=boost::python::extract<std::string>(keys[i]);
                boost::python::list tmp_values=boost::python::extract<boost::python::list>(mydict[keys[i]]);
                strVec value_str=tostrVec(tmp_values);
                convert_map[key]=value_str;
        }
        return convert_map;
}

charMap tocharmap(boost::python::dict mydict){
        charMap convert_map;
        boost::python::list keys = mydict.keys();
        for(int i = 0; i<len(keys);++i){
                std::string key=boost::python::extract<std::string>(keys[i]);
                int tmp_value=boost::python::extract<int>(mydict[keys[i]]);
                char value_str=(tmp_value==0)?'0':'1';
                convert_map[key]=value_str;
        }
        return convert_map;
}

alleleMap toallelemap(boost::python::dict mydict){
        alleleMap convert_map;
        boost::python::list keys = mydict.keys();
        for(int i = 0; i<len(keys);++i){
                std::string key=boost::python::extract<std::string>(keys[i]);
                boost::python::list tmp_values=boost::python::extract<boost::python::list>(mydict[keys[i]]);
                intVec value_str=tointVec(tmp_values);
                convert_map[key]=value_str;
        }
        return convert_map;
}

boost::python::dict todict(strdoubleMap cmap){
        boost::python::dict return_dic;
        for(auto it=cmap.begin();it!=cmap.end();it++){
                return_dic[it->first]=it->second;
        }
        return return_dic;
}

boost::python::list tolist(strVec clist){
	boost::python::list return_list;
	for (auto it=clist.begin(); it!=clist.end(); it++){
		return_list.append(*it);
	}
	return return_list;
}

boost::python::tuple apply(boost::python::list pymaf, boost::python::list invnf, boost::python::list founder, boost::python::dict parents, boost::python::dict mates, boost::python::dict offspring, boost::python::dict sex, boost::python::list sorted_miss_persons, boost::python::dict alleles, int simple){
        strVec cinvnf=tostrVec(invnf);
        strVec cfounder=tostrVec(founder);
        strstrMap cparents=tomap(parents);  
	strstrMap cmates=tomap(mates);
	strstrMap coffspring=tomap(offspring);
	charMap csex=tocharmap(sex);
	strVec csorted_miss_persons=tostrVec(sorted_miss_persons);
        alleleMap calleles=toallelemap(alleles);      
	doubleVec maf=todoubleVec(pymaf);
	mytuple returns;
        strdoubleMap conditional_prob;
	strVec missing_persons;
        boost::python::dict pyconditional_prob;
	boost::python::list pymissing_persons;
	returns=probability(&maf,&cinvnf,&cfounder,&cparents,&cmates,&coffspring,&csex,&csorted_miss_persons,&calleles,&simple);
	conditional_prob=std::get<0>(returns);
	missing_persons=std::get<1>(returns);
        pyconditional_prob=todict(conditional_prob);
	pymissing_persons=tolist(missing_persons);
        return boost::python::make_tuple(pyconditional_prob,pymissing_persons);
}

BOOST_PYTHON_MODULE(cmissing_infer){
        using namespace boost::python;
        def("probability",probability);
        def("apply",apply);
}

