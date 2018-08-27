#include <boost/python.hpp>
#include <boost/python/dict.hpp>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <iostream>
//convert int to alphabet
using namespace std;
typedef map<std::string,long double> strMap;
typedef map<std::string,strMap> MyMap;
typedef std::vector<std::string> strVec;
typedef std::vector<long double> doubleVec;

char nth_letter(int n)
{
    assert(n >= 1 && n <= 26);
    return "abcdefghijklmnopqrstuvwxyz"[n-1];
}

strVec probable(doubleVec *argmaf,strVec *arg_offspring_GT,strVec *arg_parent_GT){
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

/*int main(){
	strMap conditional_prob;
	strVec offspring_GT={"21","21"};
	long doubleVec maf={0.9,0.1};
	strVec parent_GT={"None","None"};
	conditional_prob=probable(&maf,&offspring_GT,&parent_GT);
	for (auto it=conditional_prob.begin();it!=conditional_prob.end();it++){
		cout<<it->first<<it->second<<endl;
	}
	return 0;
}*/

strVec tostrVec(boost::python::list pylist){
	strVec clist;
	for (int i=0; i<len(pylist);i++){
		string gt=boost::python::extract<std::string>(pylist[i]);
		clist.push_back(gt);
	}
	return clist;
}

doubleVec todoubleVec(boost::python::list pylist){
	doubleVec clist;
	for (int i=0; i<len(pylist);i++){
		long double gt=boost::python::extract<double>(pylist[i]);
		clist.push_back(gt);
	}
	return clist;
}

boost::python::dict todict(strMap cmap){
	boost::python::dict return_dic;
	for(auto it=cmap.begin();it!=cmap.end();it++){
		return_dic[it->first]=it->second;
	}
	return return_dic;
}

boost::python::list topythonlist(strVec pinv){
        boost::python::list postinv;
        for (auto it=pinv.begin(); it!=pinv.end(); it++){
                postinv.append(*it);
        }
        return postinv;
}

boost::python::list apply(boost::python::list pymaf, boost::python::list pyoffgt, boost::python::list pyparent_gt){
	strVec off_GT, parent_GT;
	off_GT=tostrVec(pyoffgt);
	parent_GT=tostrVec(pyparent_gt);
	doubleVec maf=todoubleVec(pymaf);
	strVec possible_GT;
	boost::python::list pypossible_GT;
	possible_GT=probable(&maf,&off_GT,&parent_GT);
	pypossible_GT=topythonlist(possible_GT);
	return pypossible_GT;
}

BOOST_PYTHON_MODULE(cmissingparents){
	using namespace boost::python;
	def("probable",probable);
	def("apply",apply);
}
