#include <boost/python.hpp>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <bitset>
#include <math.h>
using namespace std;
typedef std::vector<int> intVec;
typedef std::vector<intVec> VecintVec;
typedef std::vector<VecintVec> VecVecintVec;
typedef std::vector<std::string> strVec;
typedef map<std::string,intVec> alleleMap;
typedef map<std::string,strVec> strMap;
typedef map<std::string,std::string> simple_strMap;
typedef map<std::string,VecintVec> vMap;
typedef map<std::string,VecVecintVec> myMap;
typedef std::tuple<std::vector<strVec>, std::vector<double>> mytuple;


strVec conv(strMap*a, strVec*b)
{
	//combine the inheritance vectors of invnf
	strVec result;
	strMap tmp_dic=*a;
	strVec sorted_key=*b;
	strVec pred;
	strVec con;
	if(sorted_key.size() > 2){
		string first_key = sorted_key.front();
		strMap::iterator it = tmp_dic.find(first_key);
		pred=it->second;
		tmp_dic.erase(it);
		sorted_key.erase(sorted_key.begin());		
		con = conv(&tmp_dic,&sorted_key);
	}
	else if(sorted_key.size() == 2){
		strMap::iterator it_first = tmp_dic.find(sorted_key.at(0));
		strMap::iterator it_second = tmp_dic.find(sorted_key.at(1));
		pred=it_first->second;
		con=it_second->second;
	}
	else if(sorted_key.size() == 1){
		strMap::iterator it_first = tmp_dic.find(sorted_key.at(0));
		return it_first->second;
	}
	for (auto iter1=pred.begin(); iter1 != pred.end(); iter1++){
		for (auto iter2=con.begin(); iter2 != con.end(); iter2++){
			result.push_back(*iter1+*iter2);
		}
	}
	return result;
}

strVec postInv(strVec*arginvnf,strVec*argfounder,strMap*argparents,alleleMap*argalleles,bool rv_flag){
	strVec invnf=*arginvnf;
	strVec founder=*argfounder;
	strMap parents=*argparents;
	alleleMap alleles=*argalleles;
	strVec postinv;
	intVec tmpgt;
	for (auto it=alleles.begin(); it!=alleles.end(); it++){tmpgt.insert(tmpgt.end(),it->second.begin(),it->second.end());}
	std::sort(tmpgt.begin(),tmpgt.end());
	tmpgt.erase(std::unique(tmpgt.begin(),tmpgt.end()),tmpgt.end());
	if (tmpgt.size()==1){
		//uninformative
		const int afnfnum=invnf.size();
		int bitlen = 2*afnfnum;
		for (int it=0; it != pow(2,bitlen); it++){
			std::string s = std::bitset< 128 >( it ).to_string();
			std::string inv=s;
			int diff = s.length()-bitlen;
			assert (diff>0);
			inv.erase(inv.begin(),inv.begin()+diff);
			postinv.push_back(inv);
		}
	}
	else{
		myMap tmpv;
		vMap v;
		vMap pos;
		int offspring=0;
		strVec branch_order;
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
                                return postinv;
                        }
			strVec f_parents = {"0","0"};
			if (parents.at(tmp_parents.at(0))==f_parents && parents.at(tmp_parents.at(1))==f_parents){
				//2nd generation
				offspring = 0;
				branch_order = {*it_nf};
				for (auto it_pos=pos.at(*it_nf).begin(); it_pos!=pos.at(*it_nf).end(); it_pos++){
					intVec inh;
					intVec tmp_pos = *it_pos;
					if (tmpv.find(*it_nf) == tmpv.end()){
						tmpv.insert(std::pair<string,VecVecintVec>(*it_nf,{{}}));
					}
					else{
						tmpv.at(*it_nf).push_back({});
					}
					
					sort(tmp_pos.begin(),tmp_pos.end());
					for (auto it_ind_pos=tmp_pos.begin(); it_ind_pos!=tmp_pos.end(); it_ind_pos++){
						inh.push_back(*it_ind_pos%2);
					}
					tmpv.at(*it_nf).back().push_back(inh);
					/*if (tmpv.find(*it_nf) != tmpv.end()){
						VecintVec tmp_inh={inh};
						tmpv.at(*it_nf).push_back(tmp_inh);
					}
					else{
						VecVecintVec tmp_inh={{inh}};
						tmpv.insert(std::pair<string,VecVecintVec>(*it_nf,tmp_inh));
					}*/
				}
				for (auto it_tv=tmpv.at(*it_nf).begin(); it_tv!=tmpv.at(*it_nf).end(); it_tv++){
					VecintVec tmp_tv=*it_tv;
					if (v.find(*it_nf) != v.end()){
						v.at(*it_nf).insert(v.at(*it_nf).end(),tmp_tv.begin(),tmp_tv.end());
					}
					else{
						v.insert(std::pair<string,VecintVec>(*it_nf,tmp_tv));
					}
				}
			}
			else{
				offspring +=1;
				branch_order.push_back(*it_nf);
				string branch_id = branch_order.front();
				string fid = tmp_parents.at(0);
				string mid = tmp_parents.at(1);
				if (parents.at(fid)!=f_parents){
					//if father important
					int order_idx=find(branch_order.begin(),branch_order.end(),fid)-branch_order.begin();
					for (auto it_pos=pos.at(*it_nf).begin(); it_pos!=pos.at(*it_nf).end(); it_pos++){
						intVec element = *it_pos;
						sort(element.begin(),element.end());
						int pa = element.at(0);
						int ma = element.at(1);
						//cout<<"For pos: "<<pa<<ma<<endl;
						if (tmpv.find(*it_nf) == tmpv.end()){
							tmpv.insert(std::pair<string,VecVecintVec>(*it_nf,{{}}));
						}
						else{
							tmpv.at(*it_nf).push_back({});
						}
						for (auto it_pos_fid=pos.at(fid).begin(); it_pos_fid!=pos.at(fid).end(); it_pos_fid++){
							int idx=it_pos_fid-pos.at(fid).begin();
							intVec tmpinv;
							if ((*it_pos_fid).at(pa)<2){tmpinv={0,ma%2};}
							else if ((*it_pos_fid).at(pa)>1){tmpinv={1,ma%2};}
						  	VecintVec uniq_tmpv;
							for (auto it_tmpv=tmpv.at(fid).at(idx).begin(); it_tmpv!=tmpv.at(fid).at(idx).end();it_tmpv++){
								intVec inh;
								inh.reserve((*it_tmpv).size()+tmpinv.size());
								inh.insert(inh.end(), (*it_tmpv).begin(), (*it_tmpv).end());
								inh.insert(inh.end(), tmpinv.begin(), tmpinv.end());
								tmpv.at(*it_nf).back().push_back(inh);
								intVec fid_v;
								for (auto t_tmpv=(*it_tmpv).end()-2; t_tmpv!=(*it_tmpv).end(); t_tmpv++)
                                                                        {fid_v.push_back(*t_tmpv);}
								if (find(uniq_tmpv.begin(), uniq_tmpv.end(), fid_v)==uniq_tmpv.end()){
									uniq_tmpv.push_back(fid_v);
								}
							}
							for (auto itt_tmpv=uniq_tmpv.begin(); itt_tmpv!=uniq_tmpv.end(); itt_tmpv++){
								for (int tv_idx=0; tv_idx!=v.at(branch_id).size(); tv_idx++){
									intVec it_tv=v.at(branch_id).at(tv_idx);
									if (it_tv.size()==2*offspring){
										intVec ances_v;
										for (auto tmp_tv=it_tv.begin()+2*order_idx; tmp_tv!=it_tv.begin()+2*order_idx+2; tmp_tv++)
											{ances_v.push_back(*tmp_tv);}
										if (ances_v==*itt_tmpv){
											intVec new_v;
											new_v.reserve(it_tv.size()+tmpv.size());
											new_v.insert(new_v.end(), it_tv.begin(), it_tv.end());
											new_v.insert(new_v.end(), tmpinv.begin(), tmpinv.end());
											v.at(branch_id).push_back(new_v);
										}
									}
								}
							}
						}
					}
					for (auto it_tv=v.at(branch_id).begin(); it_tv!=v.at(branch_id).end(); it_tv++){
						if ((*it_tv).size()<2*offspring+2){v.at(branch_id).erase(it_tv);it_tv--;} //remove redundancy
					}
				}
				else{
					//if mother important
					int order_idx=find(branch_order.begin(),branch_order.end(),mid)-branch_order.begin();
					for (auto it_pos=pos.at(*it_nf).begin(); it_pos!=pos.at(*it_nf).end(); it_pos++){
						intVec element = *it_pos;
						sort(element.begin(),element.end());
						int pa = element.at(0);
						int ma = element.at(1);
						//cout<<"For pos: "<<pa<<ma<<endl;
						if (tmpv.find(*it_nf) == tmpv.end()){
							tmpv.insert(std::pair<string,VecVecintVec>(*it_nf,{{}}));
						}
						else{
							tmpv.at(*it_nf).push_back({});
						}
						for (auto it_pos_mid=pos.at(mid).begin(); it_pos_mid!=pos.at(mid).end(); it_pos_mid++){
							int idx=it_pos_mid-pos.at(mid).begin();
							intVec tmpinv;
							if ((*it_pos_mid).at(ma-2)<2){tmpinv={pa%2,0};}
							else if ((*it_pos_mid).at(ma-2)>1){tmpinv={pa%2,1};}
							VecintVec uniq_tmpv;
							for (auto it_tmpv=tmpv.at(mid).at(idx).begin(); it_tmpv!=tmpv.at(mid).at(idx).end();it_tmpv++){
								intVec inh;
								inh.reserve((*it_tmpv).size()+tmpinv.size());
								inh.insert(inh.end(), (*it_tmpv).begin(), (*it_tmpv).end());
								inh.insert(inh.end(), tmpinv.begin(), tmpinv.end());
								tmpv.at(*it_nf).back().push_back(inh);
                                                                intVec mid_v;
                                                                for (auto t_tmpv=(*it_tmpv).end()-2; t_tmpv!=(*it_tmpv).end(); t_tmpv++)
                                                                        {mid_v.push_back(*t_tmpv);}
                                                                if (find(uniq_tmpv.begin(), uniq_tmpv.end(), mid_v)==uniq_tmpv.end()){
                                                                        uniq_tmpv.push_back(mid_v);
                                                                }
							}
							for (auto itt_tmpv=uniq_tmpv.begin(); itt_tmpv!=uniq_tmpv.end(); itt_tmpv++){
								int v_size = v.at(branch_id).size();
								for (int tv_idx=0; tv_idx!=v_size; tv_idx++){
									intVec it_tv = v.at(branch_id).at(tv_idx);
									if (it_tv.size()==2*offspring){
										intVec ances_v;
										for (auto tmp_tv=it_tv.begin()+2*order_idx; tmp_tv!=it_tv.begin()+2*order_idx+2;tmp_tv++)
											{ances_v.push_back(*tmp_tv);}
										if (ances_v==*itt_tmpv){
											intVec new_v;
											new_v.reserve(it_tv.size()+tmpv.size());
											new_v.insert(new_v.end(), it_tv.begin(), it_tv.end());
											new_v.insert(new_v.end(), tmpinv.begin(), tmpinv.end());
											v.at(branch_id).push_back(new_v);
										}
									}
								}
							}
						}
					}
					for (auto it_tv=v.at(branch_id).begin(); it_tv!=v.at(branch_id).end(); it_tv++){
						if ((*it_tv).size()<2*offspring+2){v.at(branch_id).erase(it_tv);it_tv--;} //remove redundancy
					}
				}

			}
			
		}
		strMap trans_v;
		strVec sorted_key;
		for (auto it_v=v.begin(); it_v!=v.end(); it_v++){
			strVec tmp_v;
			sorted_key.push_back(it_v->first);
			for (auto it_vecint=it_v->second.begin(); it_vecint!=it_v->second.end(); it_vecint++){
				string str_v;
				for (auto it_int=(*it_vecint).begin(); it_int!=(*it_vecint).end(); it_int++){
					str_v += std::to_string(*it_int);
				}
				tmp_v.push_back(str_v);
			}
			trans_v.insert(std::pair<string,strVec>(it_v->first,tmp_v));
		}	
		std::sort(sorted_key.begin(),sorted_key.end(),[invnf](const string &a, const string &b)
		{
			return (find(invnf.begin(),invnf.end(),a) < find(invnf.begin(),invnf.end(),b));
		});
		postinv=conv(&trans_v,&sorted_key);
		if (rv_flag){
			intVec rv_idx;
			simple_strMap rv_invs;
			intVec one={1,1};
			string first_pv=postinv.at(0);
			for (auto it_invnf=invnf.begin(); it_invnf!=invnf.end(); it_invnf++){
				int invnf_idx=it_invnf-invnf.begin();
				if (alleles.at(*it_invnf)!=one){
					rv_idx.push_back(invnf_idx);
				}
			}	
			for (auto it_pv=postinv.begin(); it_pv!=postinv.end(); it_pv++){
				string rv_inv;
				for (auto it_rv_idx=rv_idx.begin(); it_rv_idx!=rv_idx.end(); it_rv_idx++){
					rv_inv+=(*it_pv).substr(2*(*it_rv_idx),2);
				}
				if(rv_invs.find(rv_inv)==rv_invs.end()){
					rv_invs.insert(std::pair<string,string>(rv_inv,*it_pv));
				}
				else{
					*it_pv=rv_invs.at(rv_inv);
				}
			}
		}
	}
	return postinv;	
}

strMap assign_allele(strVec *arginvnf, strVec *arginvnf_sorted, strVec *arg_founder, string *arg_v, strMap *arg_parents){
	strVec invnf = *arginvnf;
	strVec invnf_sorted = *arginvnf_sorted;
	strVec founder = *arg_founder;
	string v=*arg_v;
	strMap parents= *arg_parents;
	strMap fam_allele;
	for (auto it=invnf_sorted.begin(); it!=invnf_sorted.end(); it++){
		int idx = find(invnf.begin(),invnf.end(),*it)-invnf.begin();
		int mv1=v.at(2*idx)-'0';
		int mv2=v.at(2*idx+1)-'0';
		string fid = parents.at(*it).at(0);
		string mid = parents.at(*it).at(1);
		strVec allele;
		if (fam_allele.find(fid)!=fam_allele.end()){
			allele.push_back(fam_allele.at(fid).at(mv1));
		}
		else {allele.push_back(fid+v.at(2*idx));}
		if (fam_allele.find(mid)!=fam_allele.end()){
			allele.push_back(fam_allele.at(mid).at(mv2));
		}
		else {allele.push_back(mid+v.at(2*idx+1));}
		fam_allele.insert(std::pair<string,strVec>(*it,allele));
	}
	for (auto it=founder.begin(); it!=founder.end(); it++){
		strVec allele={*it+"0",*it+"1"};
		fam_allele.insert(std::pair<string,strVec>(*it,allele));
	}
	return fam_allele;	
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

strMap tomap(boost::python::dict mydict){
        strMap convert_map;
        boost::python::list keys = mydict.keys();
        for(int i = 0; i<len(keys);++i){
                std::string key=boost::python::extract<std::string>(keys[i]);
		boost::python::list tmp_values=boost::python::extract<boost::python::list>(mydict[keys[i]]);
		strVec value_str=tostrVec(tmp_values);
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

boost::python::list topythonlist(strVec pinv){
	boost::python::list postinv;
	for (auto it=pinv.begin(); it!=pinv.end(); it++){
		postinv.append(*it);
	}
	return postinv;
}

boost::python::dict todict(strMap allele){
	boost::python::dict tallele;
	for (auto it=allele.begin(); it!=allele.end(); it++){
		boost::python::list a;
		for (auto itt=it->second.begin(); itt!=it->second.end(); itt++){
			a.append(*itt);
		}
		tallele[it->first]=a;
	}
	return tallele;
}

boost::python::list apply(boost::python::list invnf, boost::python::list founder, boost::python::dict parents, boost::python::dict alleles, bool rv_flag){
	strVec cinvnf=tostrVec(invnf);
	strVec cfounder=tostrVec(founder);
	strMap cparents=tomap(parents);
	alleleMap calleles=toallelemap(alleles);
	strVec cpostinv;
	cpostinv=postInv(&cinvnf,&cfounder,&cparents,&calleles,rv_flag);
	boost::python::list py_postinv;
	py_postinv=topythonlist(cpostinv);
	return py_postinv;
}

boost::python::dict apply_assign(boost::python::list invnf, boost::python::list invnf_sorted, boost::python::list founder, string v, boost::python::dict parents){
	strVec cinvnf=tostrVec(invnf);
	strVec cinvnf_sorted=tostrVec(invnf_sorted);
	strVec cfounder=tostrVec(founder);
	strMap cparents=tomap(parents);
	strMap fam_allele;
	fam_allele=assign_allele(&cinvnf,&cinvnf_sorted,&cfounder,&v, &cparents);
	boost::python::dict py_fam_allele=todict(fam_allele);
	return py_fam_allele;
}


BOOST_PYTHON_MODULE(cpostInv){
        using namespace boost::python;
        def("apply",apply);
        def("apply_assign",apply_assign);
}

/*int main(){
	strVec invnf = {"3","5","6","9","11","12"};
	strVec founder = {"1","2","4","10"};
	strMap parents = {{"1",{"0","0"}},{"2",{"0","0"}},{"4",{"0","0"}},{"10",{"0","0"}},{"3",{"1","2"}},{"5",{"4","3"}},{"6",{"4","3"}},{"9",{"1","2"}},{"11",{"9","10"}},{"12",{"9","10"}}};
	alleleMap alleles = {{"1",{1,2}},{"2",{1,2}},{"3",{2,2}},{"4",{1,1}},{"5",{1,2}},{"6",{1,2}},{"9",{2,2}},{"10",{1,1}},{"11",{1,2}},{"12",{1,2}}};
	strVec postinv=postInv(&invnf,&founder,&parents,&alleles);
	//for (auto it=postinv.begin(); it!=postinv.end(); it++){
	//	cout<<*it<<endl;
	//}
	//strMap a = {{"3",{"1001","1000","1010","1011"}},
	//	    {"9",{"0101","0100","0110","0111"}}};
	//strVec sorted_key = {"3","9"};
	//strVec postinv = conv(&a,&sorted_key);
	//for (auto it=postinv.begin(); it != postinv.end(); it++){
	//	std::cout<<*it<<" "<<std::endl;
	//}
}*/

