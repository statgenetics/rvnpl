#include <boost/python.hpp>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include <stdlib.h>
#include <sstream>
#include <algorithm>
#include <iostream>
using namespace std;
typedef map<std::string,double> strMap;
typedef map<std::string,strMap> MyMap;
typedef std::vector<std::string> strVec;
typedef std::tuple<std::vector<strVec>, std::vector<double>> mytuple;
mytuple conv(MyMap*a, float min_prob)
{
	MyMap comb;
	MyMap tmp_dic;
	strMap first_dict;
	tmp_dic = *a;
	strVec key_list;
	std::vector<strVec> return_gt;
	std::vector<double> return_prob;
	for(auto iter=tmp_dic.begin();iter != tmp_dic.end();iter++){
		key_list.push_back(iter->first);
	}
	std::sort(key_list.begin(),key_list.end(),std::greater<std::string>());
	string first_key = key_list.front();
	key_list.pop_back();
	MyMap::iterator it=tmp_dic.find(first_key);
	for (auto inner_it=it->second.begin();inner_it!=it->second.end();inner_it++){
		first_dict[inner_it->first]=inner_it->second;
	}
	tmp_dic.erase(it);
	if(tmp_dic.empty()){
		//last one
		for (auto gene_it=first_dict.begin();gene_it!=first_dict.end();gene_it++){
			if (gene_it->second > min_prob){
				strVec new_gene_config;
				new_gene_config.push_back(gene_it->first);
				return_gt.push_back(new_gene_config);
				return_prob.push_back(gene_it->second);
			}
		}
	}
	else{
		std::vector<strVec> pre_gt;
		std::vector<double> pre_prob;
		mytuple pre_tuple=conv(&tmp_dic, min_prob);
		pre_gt=std::get<0>(pre_tuple);
		pre_prob=std::get<1>(pre_tuple);
		for(auto pre_it=pre_gt.begin();pre_it!=pre_gt.end();pre_it++){
			for(auto first_it=first_dict.begin();first_it!=first_dict.end();first_it++){
				double comb_prob;
				int pos = pre_it-pre_gt.begin();
				comb_prob=pre_prob.at(pos)*first_it->second;
				if(comb_prob>min_prob){
					strVec combined_gt=pre_gt.at(pos);
					combined_gt.push_back(first_it->first);
					return_gt.push_back(combined_gt);
					return_prob.push_back(comb_prob);
				}
			}
		}
	}
	double sum_prob=0;
	for(auto prob_it=return_prob.begin();prob_it!=return_prob.end();prob_it++){
		sum_prob+=*prob_it;
	}
	//cout<<"sum_prob"<<sum_prob<<endl;
	/*if(sum_prob!=1){
		for(auto prob_it=return_prob.begin();prob_it!=return_prob.end();prob_it++){
			int pos=prob_it-return_prob.begin();
			return_prob.at(pos)=*prob_it/sum_prob;
		}
	}*/
	mytuple returns = std::make_tuple(return_gt,return_prob);
	return returns;
}

int main(){
	strMap off1 = {
		{"12",0.2},{"11",0.7},{"13",0.1}
	};
	strMap off2 = {
		{"12",0.3},{"11",0.6},{"23",0.1}
	};
	strMap off3 = {
		{"12",0.3},{"11",0.6},{"22",0.1}
	};
	MyMap dic = {
		{"off1",off1},{"off2",off2},{"off3",off3}
	};
	for(auto it=dic.begin();it!=dic.end();it++){
		for(auto inner_it=it->second.begin();inner_it!=it->second.end();inner_it++){
		}
	}
	mytuple combined=conv(&dic, 1E-5);
	std::vector<strVec> combined_gt=std::get<0>(combined);
	std::vector<double> comb_prob=std::get<1>(combined);
	for(auto gt_it=combined_gt.begin();gt_it!=combined_gt.end();gt_it++){
		int pos = gt_it-combined_gt.begin();
		cout<<comb_prob.at(pos)<<endl;
		for(auto inner_it=gt_it->begin();inner_it!=gt_it->end();inner_it++){
			cout<<*inner_it;
		}
		cout<<endl;
	}
	return 0;
}

strMap tomap(boost::python::dict mydict){
	strMap convert_map;
	//boost::python::extract<boost::python::dict> cppdict_extract(pre_dic);
	boost::python::list keys = mydict.keys();
	for(int i = 0; i<len(keys);++i){
		std::string gt=boost::python::extract<std::string>(keys[i]);
		double prob=boost::python::extract<double>(mydict[keys[i]]);
		convert_map[gt]=prob;
	}
	
	return convert_map;
}

boost::python::list togtlist(std::vector<strVec> *gt){
  	std::vector<strVec> tmpgt = *gt;
	boost::python::list gt_list;
	for (auto iter=tmpgt.begin();iter != tmpgt.end();iter++){
		boost::python::list tmp_list;
		for (auto inner_it=iter->begin();inner_it!=iter->end();inner_it++){
			tmp_list.append(*inner_it);
		}
		gt_list.append(tmp_list);
	}
	return gt_list;
}

boost::python::list toproblist(std::vector<double> *prob){
  	std::vector<double> tmpprob = *prob;
	boost::python::list prob_list;
	for (auto iter=tmpprob.begin();iter != tmpprob.end();iter++){
		prob_list.append(*iter);
	}
	return prob_list;
}

boost::python::tuple applyconv(boost::python::dict dic, float min_prob = 1E-5){
	MyMap map1;
	boost::python::dict mydict;
	mytuple returns;
	std::vector<strVec> return_gt;
	std::vector<double> return_prob;
	boost::python::list py_return_gt;
	boost::python::list py_return_prob;
	boost::python::list keys = dic.keys();
	for(int i=0; i<len(keys);i++){
		std::string dkey=boost::python::extract<std::string>(keys[i]);
		strMap tmp_dic;
		boost::python::extract<boost::python::dict> inner_dic1(dic[dkey]);
		if(inner_dic1.check()){
			tmp_dic=tomap(inner_dic1);
			map1[dkey]=tmp_dic;
		}
	}
	returns = conv(&map1, min_prob);
	return_gt = std::get<0>(returns);
	return_prob=std::get<1>(returns);
	py_return_gt = togtlist(&return_gt);
	py_return_prob = toproblist(&return_prob);
	return boost::python::make_tuple(py_return_gt,py_return_prob);
}

BOOST_PYTHON_MODULE(cconv2){
	using namespace boost::python;
	def("conv",conv);
	def("applyconv",applyconv);
}
