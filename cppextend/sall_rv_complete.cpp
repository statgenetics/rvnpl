//
//  sall.cpp
//  Percy
//
//  Created by Linhai Zhao on 4/22/15.
//  Copyright (c) 2015 Linhai Zhao. All rights reserved.
//

#include <boost/python.hpp>
#include <math.h>
#include <map>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
typedef std::vector<int> IntVector;
typedef std::vector<std::string> StringVector;
typedef std::vector<StringVector> VecStringVec;
typedef map<int, StringVector> MyMap;


/*int main(){
	double s_all;
    StringVector fallele;
    MyMap fam_allele;
    IntVector affnf;
    affnf={2,3};
    fallele={"gp0","gp1","gm0","gm1"};
    fam_allele[2]={"gp0","gm0"};
    fam_allele[3]={"gp1","gm0"};
    fam_allele[4]={"gp0","gm1"};
    s_all=sall(&affnf,&fam_allele, &fallele);
    return 0;
}*/

int factorial(int n){
    int fac=1;
    if (n==0) {
        return 1;
    }
    if (n<0) {
        return -99;
    }
    while (n>0) {
        fac = fac*n;
        n--;
    }
    return fac;
}

double sall(IntVector*aff, MyMap*fam_a,StringVector*fa, IntVector*fm,bool rv_flag){
    IntVector tmp_aff = *aff;
    int affnum = tmp_aff.size();
    double sall;
    MyMap fam_allele = *fam_a;
    StringVector founderallele=*fa;
    IntVector foundermarker=*fm;
    int total = pow(2.0,affnum);
    int total_count=0;
    for (int i=0; i<total; i++) {
        //for each unique selection of alleles
        //to binary;
        string bin;
        int number = i;
        char holder=' ';
        if (number==0) {
            bin='0';
        }
        while(number !=0){
            holder=number%2+'0';
            bin=holder+bin;
            number /=2;
        }
        //cout<<string((affnum-bin.size()),'0')<<endl;
        bin=string((affnum-bin.size()),'0')+bin;
        StringVector pick(affnum);
        //convert binary digits to founder alleles
        for (int j=0; j<bin.size(); j++) {
            int binj = bin[j]-'0'; //convert character to integer
            pick[j]=fam_allele[tmp_aff[j]][binj];
        }
        //h[i]=pick;
        int tmp_count = 1;
        //count the occurrences of founder alleles in picked alleles
        for (StringVector::iterator a=founderallele.begin(); a != founderallele.end(); ++a) {
            int occurrence=0;// = std::count(pick.begin(), pick.end(), *a);
	    if (rv_flag && foundermarker.at(a-founderallele.begin()) != 1 || !rv_flag){
			for (StringVector::iterator pit=pick.begin(); pit != pick.end();++pit){
				if (*pit == *a){occurrence += 1;}
				}
			}
            tmp_count *= factorial(occurrence);
        }
        total_count += tmp_count;
    }
    sall=double(total_count)/double(total);
	return sall;
}

IntVector toIntVec(boost::python::list aff){
	IntVector t_aff;
	for(int i =0;i<len(aff);++i){
		int tmp=boost::python::extract<int>(aff[i]);
		t_aff.push_back(tmp);
	}
	return t_aff;
}

StringVector toStrVec(boost::python::list fa){
	StringVector t_fa;
	for(int i = 0; i<len(fa);++i){
		string tmp=boost::python::extract<string>(fa[i]);
		t_fa.push_back(tmp);
	}
	return t_fa;
}

MyMap tomap(boost::python::dict nf){
	MyMap t_nf;
	boost::python::list keys = nf.keys();
	for(int i=0; i<len(keys);++i){
		int extracted_key=boost::python::extract<int>(keys[i]);
		StringVector extracted_val;
		extracted_val.push_back(boost::python::extract<string>(nf[extracted_key][0]));
		extracted_val.push_back(boost::python::extract<string>(nf[extracted_key][1]));
		t_nf[extracted_key]=extracted_val;
	}
	return t_nf;
}

double apply(boost::python::list aff, boost::python::dict fam_a, boost::python::list fa, boost::python::list fm, bool rv_flag){
	IntVector caff;
	StringVector fallele;
	IntVector fmarker;
	MyMap cfam_a;
	double s_all;
	caff = toIntVec(aff);
	fallele = toStrVec(fa);
	fmarker = toIntVec(fm);
	cfam_a = tomap(fam_a);
	s_all = sall(&caff,&cfam_a,&fallele,&fmarker,rv_flag);
	return s_all;
}

BOOST_PYTHON_MODULE(sall_cpp){
	using namespace boost::python;
	def("sall",sall);
	def("apply",apply);
}
