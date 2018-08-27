//ibd_m.cpp equivalent cpp script of ibd_m.py
#include <boost/python.hpp>
#include<iostream>
#include<vector>
#include<algorithm>
#include<map>
#include<string>
using namespace std;
typedef std::vector<float> VecFloat;
typedef std::vector<string> VecString;
typedef std::vector<int> VecInt;
typedef std::vector<int> VecInt;
typedef std::vector<VecInt> VecVecInt;
VecFloat sib_ibd(VecString &geno);
float prob_ibd(VecString &geno, string &allele);
float cousin_ibd(VecString &geno);
float un_ibd(VecString &geno); 

VecFloat sib_ibd(VecString &geno)
{
	VecFloat ibd (3,0.0);
	VecVecInt inherit;
	std::vector<VecVecInt> compat_inherit;
	for (VecString::iterator it=geno.begin()+4; it<geno.end(); it++)
	{
		//for each allele in offspring
		string sib_allele = *it;
		VecInt tmp_inheirt;
		VecString::iterator pos=geno.begin()-1;
		while (pos != geno.begin()+4)
		{
			pos = std::find(pos+1,geno.begin()+4,sib_allele);
			if ((pos-geno.begin()) !=4){tmp_inheirt.push_back(pos-geno.begin());}
		}
		inherit.push_back(tmp_inheirt);
	}
	for (int indiv=0; indiv<2; indiv++)
	{
		//for each offspring
		VecVecInt sib_inherit;
		for (VecInt::iterator tmp_it=inherit[2*indiv].begin();tmp_it<inherit[2*indiv].end();tmp_it++)
		{
			for (VecInt::iterator tmp2_it=inherit[2*indiv+1].begin();tmp2_it<inherit[2*indiv+1].end();tmp2_it++)
			{
				if (*tmp_it<2 && *tmp2_it>1 || *tmp_it>1 && *tmp2_it<2)
				{
					VecInt temp {*tmp_it,*tmp2_it};
					sib_inherit.push_back(temp);
				}
			}
		}
		compat_inherit.push_back(sib_inherit);
	}
	float count[3] = {0.0,0.0,0.0};
	int count_sum = compat_inherit[0].size()*compat_inherit[1].size();
	for (VecVecInt::iterator tmp1_it=compat_inherit[0].begin(); tmp1_it<compat_inherit[0].end(); tmp1_it++)
	{
		//count IBD value
		for (VecVecInt::iterator tmp2_it=compat_inherit[1].begin(); tmp2_it<compat_inherit[1].end(); tmp2_it++)
		{
			int ibd_count=0;
			VecInt tmp_common;
			for (VecInt::iterator i = tmp1_it->begin(); i != tmp1_it->begin()+2;i++)
			{
				if (std::find(tmp2_it->begin(),tmp2_it->begin()+2,*i) != tmp2_it->begin()+2 && geno.at(*i) != "1")
				{
					tmp_common.push_back(*i);//shared allele in [0,1,2,3]
				}
			}
			ibd_count = tmp_common.size();
			count[ibd_count]++;
		}
	}
	for (int i=0; i<3; i++){ibd[i]=count[i]/count_sum;}
	return ibd;	
}

float prob_ibd(VecString &geno, string &allele)
{
	VecInt pos;   //position of the given allele in parental genotypes
	VecVecInt alt_pos;
	float prob=0.0;
	VecString alt {"-1","-1"};
	for (auto iter=geno.begin(); iter != geno.begin()+4; iter++)
	{
		if (*iter == allele){pos.push_back((iter-geno.begin()));}
	}
	for (auto iter=geno.begin()+4; iter != geno.end(); iter++)
	{
		//determine the alternative allele
		int tmp_pos = iter-(geno.begin()+4);
		if (*iter != allele)
		{
			if (tmp_pos<2){alt[0]=*iter;}
			else {alt[1]=*iter;}
		}
		else if (tmp_pos%2==1 && alt[tmp_pos/2]=="-1"){alt[tmp_pos/2]=*iter;}
	}
	VecInt altpos1, altpos2;
	for (auto iter=geno.begin(); iter != geno.begin()+4; iter++)
	{
		if (*iter == alt[0]){altpos1.push_back((iter-geno.begin()));}
		if (*iter == alt[1]){altpos2.push_back((iter-geno.begin()));}
	}
	alt_pos.push_back(altpos1);
	alt_pos.push_back(altpos2);
	vector<map<int,float>> freq_map;
	for (int i=0;i<2;i++)
	{
		map<int,float> tmp_occurrence;
		VecInt tmp_p_compat;
		for (auto p : pos)
		{
			for (auto q : alt_pos[i])
			{
				if (p<2 && q>1 || p>1 && q<2){tmp_p_compat.push_back(p);}
			}
		}
		for (auto tmp : tmp_p_compat)
		{
			float count = (float) 1/tmp_p_compat.size();
			if(tmp_occurrence.find(tmp)==tmp_occurrence.end())
			{
				//not in key
				tmp_occurrence[tmp]=count;
			}
			else {tmp_occurrence[tmp]+=count;}
		}
		freq_map.push_back(tmp_occurrence);
	}
	for (auto &it1 : freq_map[0])
	{
		for (auto &it2 : freq_map[1])
		{
			if (it2.first==it1.first){prob+=it1.second*it2.second;}
		}
	}
	return prob;
}

float cousin_ibd(VecString &geno)
{
	//calculate IBD between cousins
	//geno = GT for [grandparents, fam1, fam2] important parent put in the first place of each family
	VecString ag = {geno[8],geno[9]};       //GT for 2 cousins
	VecString bg = {geno[14],geno[15]};
	VecString shared_allele;
	float total_p = 0.0;
	for (auto ait : ag)
	{
		if (std::find(bg.begin(),bg.end(),ait) != bg.end() && ait != "1"){shared_allele.push_back(ait);} //possible shared RV
	}
	if (shared_allele.size()==2 && shared_allele[0]==shared_allele[1]){shared_allele.pop_back();}
	for (auto a : shared_allele)
	{	
		int rep=0;
		for (auto tmp : ag)
		{
			for (auto tmp1 : bg)
			{
				if (tmp==a && tmp1==a){rep++;}
			}
		}
		VecFloat f_inherit;
		int flag=0;
		for (int fid=0; fid<2; fid++)
		{
			VecVecInt possible_inherit;
			VecString pg = {geno[fid*6+4],geno[fid*6+5],geno[fid*6+6],geno[fid*6+7]};     //parental genotypes
			string alt="-1";
			VecInt pos, pos_alt;
			if (fid==0)
			{
				for (auto ait=ag.begin(); ait != ag.end(); ait++)
				{
					if (*ait != a)
					{
						alt=*ait;
						break;
					}
					else if (ait-ag.begin()==1)
					{
						alt=*ait;
					}
				}
			}
			if (fid==1)
			{
				for (auto bit=bg.begin(); bit != bg.end(); bit++)
				{
					if (*bit != a)
					{
						alt=*bit;
						break;
					}
					else if (bit-bg.begin()==1)
					{
						alt=*bit;
					}
				}
			}
			for (auto it = pg.begin(); it !=pg.end(); it++)
			{
				if (*it == a) {pos.push_back((it-pg.begin()));}
				if (*it == alt) {pos_alt.push_back((it-pg.begin()));}
			}
			for (auto p : pos)
			{
				for (auto q : pos_alt)
				{
					if (p<2 && q>1 || p>1 && q<2)
					{
						VecInt tmp_pos={p,q};
						//cout<<p<<" "<<q<<endl;
						possible_inherit.push_back(tmp_pos);
					}
				}
			}
			int count=0;
			for (auto inherit : possible_inherit)
			{
				if (inherit[0]<2){count++;}
			}
			float f_in;
			if (possible_inherit.size()>0)
			{
				f_in = (float)count/possible_inherit.size();
			}
			else {f_in=0;}
			f_inherit.push_back(f_in);
			if (f_in>0){flag++;}     //family having a possibility to inherit the particular allele
			if (flag==2)
			{
				//get the GT for upper nuclear family
				VecString up_gt={geno[0],geno[1],geno[2],geno[3],geno[4],geno[5],geno[10],geno[11]};  
				float prob = prob_ibd(up_gt,a);
				total_p += (f_inherit[0]*f_inherit[1]*prob)*rep;
			}
			else {total_p+=0;}
		}
	}
	return total_p;
}

float un_ibd(VecString &geno)
{
	//calculate IBD between Uncle-Nephew pair
	//geno = GT for [grandparents uncle father mother kid(nephew)]
	VecString Uncle_g = {geno[4],geno[5]};
	VecString Nephew_g = {geno[10],geno[11]};
	VecString pg = {geno[6],geno[7],geno[8],geno[9]};
	VecString gt = {geno[0],geno[1],geno[2],geno[3],geno[4],geno[5],geno[6],geno[7]};
	VecString shared_allele;
	float total_p=0.0;
	for (auto g : Uncle_g)
	{
		if (std::find(Nephew_g.begin(),Nephew_g.end(),g) != Nephew_g.end() && g != "1")
		{
			shared_allele.push_back(g);
		}
	}
	if (shared_allele.size()==2 && shared_allele[0]==shared_allele[1]){shared_allele.pop_back();}
	for (auto ref : shared_allele)
	{
		VecString tmp_alleles;
		for (auto b : Nephew_g)
		{
			if (b==ref){tmp_alleles.push_back(b);}
		}
		for (auto allele : tmp_alleles)
		{
			float f=0.0;
			int m=0;
			int n=0;
			string alt="-1";
			for (auto bit=Nephew_g.begin(); bit != Nephew_g.end(); bit++)
			{
				if (*bit != ref)
				{
					alt=*bit;
					break;
				}
				else if (bit-Nephew_g.begin()==1)
				{
					alt=*bit;
				}
			}
			VecInt pos1,pos2;
			for (auto tmp=pg.begin();tmp!=pg.end();tmp++)
			{
				if (*tmp==ref){pos1.push_back((tmp-pg.begin()));}
				if (*tmp==alt){pos2.push_back((tmp-pg.begin()));}
			}
			for (auto p : pos1)
			{
				for (auto q : pos2)
				{
					if (p<2 && q>1){m++;}
					else if (p>1 && q<2){n++;}
				}
			}
			if ((m+n)>0){f=(float)m/(m+n);}
			else{f=0;}
			float prob=0;
			if (f!=0)
			{
				prob=prob_ibd(gt,ref);
			}
			if (Uncle_g[0]==Uncle_g[1] && Uncle_g[0]==ref)
			{
				prob *=2;
			}
			total_p += f*prob;
		}
	}
	return total_p;
}
boost::python::list sib_apply(boost::python::list geno_py)
{
	VecString geno;
	VecFloat ibd;
	boost::python::list ibd_py;
	for (int i=0; i<len(geno_py); i++)
	{
		boost::python::extract<int> extracted_geno(geno_py[i]);
		string str_geno = std::to_string(extracted_geno); 
		geno.push_back(str_geno);
	} 
	ibd=sib_ibd(geno);
	for (auto iter : ibd)
	{
		ibd_py.append(iter);
	}
	return ibd_py;
}

float prob_ibd_apply(boost::python::list geno_py,int allele)
{
	VecString geno;
	float prob;
	string str_allele = std::to_string(allele);
	for (int i=0; i<len(geno_py); i++)
	{
		boost::python::extract<int> extracted_geno(geno_py[i]);
		string str_geno = std::to_string(extracted_geno);
		geno.push_back(str_geno);
	} 
	prob = prob_ibd(geno, str_allele);
	return prob;
}

float cousin_apply(boost::python::list geno_py)
{
	VecString geno;
	float prob;
	for (int i=0; i<len(geno_py); i++)
	{
		boost::python::extract<int> extracted_geno(geno_py[i]);
		string str_geno = std::to_string(extracted_geno);
		geno.push_back(str_geno);
	} 
	prob = cousin_ibd(geno);
	return prob;
}

float un_apply(boost::python::list geno_py)
{
	VecString geno;
	float prob;
	for (int i=0; i<len(geno_py); i++)
	{
		boost::python::extract<int> extracted_geno(geno_py[i]);
		string str_geno = std::to_string(extracted_geno);
		geno.push_back(str_geno);
	} 
	prob = un_ibd(geno);
	return prob;
}
BOOST_PYTHON_MODULE(ibd_rv_cpp){
	using namespace boost::python;
	def("sib_ibd",sib_ibd);
	def("sib_apply",sib_apply);
	def("prob_ibd_apply",prob_ibd_apply);
	def("cousin_apply",cousin_apply);
	def("un_apply",un_apply);
}

/*int main()
{
	VecString geno {"1","1","2","1","1","1","2","1","1","1","1","1"};
	float prob=un_ibd(geno);
	cout<<"IBD:"<<prob<<endl;
	return 0;
}
int main(int argc, char* argv[])
{
	VecString geno;
	VecFloat ibd;
	for (char* i=*(argv+1); *i != '\0';i++)
	{
		geno.push_back(string(1,*i));
	}
	ibd=sib_ibd(geno);
}*/
