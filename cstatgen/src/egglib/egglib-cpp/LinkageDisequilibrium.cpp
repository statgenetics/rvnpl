/*
    Copyright 2009,2011 Stephane De Mita, Mathieu Siol
    
    This file is part of the EggLib library.

    EggLib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    EggLib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with EggLib.  If not, see <http://www.gnu.org/licenses/>.
  
*/


#include "LinkageDisequilibrium.hpp"
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <vector>

namespace egglib {

    // helpers

    void LinkageDisequilibrium::init() {
        BaseDiversity::init();
        _n = 0;
        _d = NULL;
        _D = NULL;
        _Dp = NULL;
        _r = NULL;
        _r2 = NULL;
        _site1 = NULL;
        _site2 = NULL;
    }

    void LinkageDisequilibrium::reset() {
        _d  = (int*)    realloc(_d , sizeof(int)    * _n);
        _D  = (double*) realloc(_D , sizeof(double) * _n);
        _Dp = (double*) realloc(_Dp, sizeof(double) * _n);
        _r  = (double*) realloc(_r , sizeof(double) * _n);
        _r2 = (double*) realloc(_r2, sizeof(double) * _n);
        _site1 = (unsigned int*) realloc(_site1, sizeof(unsigned int) * _n);
        _site2 = (unsigned int*) realloc(_site2, sizeof(unsigned int) * _n);

        if (!(_d && _D && _Dp && _r && _r2 && _site1 && _site2)) {
            clear();
            init();
            throw EggMemoryError();
        }
    }


    void LinkageDisequilibrium::clear() {
        BaseDiversity::clear();
        if (_d)  free(_d);
        if (_D)  free(_D);
        if (_Dp) free(_Dp);
        if (_r)  free(_r);
        if (_r2) free(_r2);
        if (_site1) free(_site1);
        if (_site2) free(_site2);
    }

    // constructors

    LinkageDisequilibrium::LinkageDisequilibrium() {
        init();
    }


    LinkageDisequilibrium::~LinkageDisequilibrium() { 
        clear();
    }


    // main method
    
    void LinkageDisequilibrium::load(CharMatrix& data,
                double minimumExploitableData,  unsigned int ignoreFrequency,
                std::string characterMapping) {
                    
        clear();
        init();
        importSites(data, false, minimumExploitableData, ignoreFrequency,
                                        characterMapping, false, false);
                    
        for (unsigned int i=0; i<v_S; i++) {
            if (v_sites[i]->numberOfAlleles()!=2) continue;
            for (unsigned int j=i+1; j<v_S; j++) {
                if (v_sites[j]->numberOfAlleles()!=2) continue;
                add(data, i, j);
            }
        }
    }


    // iterative filler

    void LinkageDisequilibrium::add(CharMatrix& data, unsigned int position1, unsigned int position2) {

        // initializes stuff
        unsigned int A  = 0;
        unsigned int B  = 0;
        unsigned int AB = 0;
        unsigned int c  = 0;

        char cA = v_sites[position1]->allele(0);
        char cB = v_sites[position2]->allele(0);
        
        // counts all possible alleles

        for (unsigned int i=0; i<v_ns; i++) {
            if (data.populationLabel(i)==999) continue;

            char char1 = data.character(i, v_sitePositions[position1]);
            if (!(p_characterMapping.find(char1) < p_pos_sep_mapping)) continue;
            
            char char2 = data.character(i, v_sitePositions[position2]);
            if (!(p_characterMapping.find(char2) < p_pos_sep_mapping)) continue;

            c++;
            if (char1==cA) {
                A++;
                if (char2==cB) {
                    B++;
                    AB++;
                }
            }
            else if (char2==cB) B++;
        }

        if (!c) return;

        // computes

        double pA  = (float)A/c;
        double pB  = (float)B/c;
        double pAB = (float)AB/c;
        double D = pAB - pA*pB;

        double Dmin = max( pA*pB, (1-pA)*(1-pB) ); 
        double Dmax = min( pA*(1-pB), (1-pA)*pB );
        double Dp = D/Dmax; // if D>=0
        if (D<0) Dp = D/Dmin;
        double r = D / (sqrt(pA*(1-pA)*pB*(1-pB)));
        double r2 = r*r;

        // adding the entry
        _n++;
        reset();

        _site1[_n-1] = data.sitePosition( v_sitePositions[position1] );
        _site2[_n-1] = data.sitePosition( v_sitePositions[position2] );
        _d[_n-1] = _site2[_n-1] - _site1[_n-1];
        _D[_n-1]  = D;
        _Dp[_n-1] = Dp;
        _r[_n-1]  = r;
        _r2[_n-1] = r2;

    }


    unsigned int LinkageDisequilibrium::numberOfPairs() const {
        return _n;
    }

    int LinkageDisequilibrium::d(unsigned int pair_index) {
        check(pair_index);
        return _d[pair_index];
    }

    double LinkageDisequilibrium::D(unsigned int pair_index) {
        check(pair_index);
        return _D[pair_index];
    }

    double LinkageDisequilibrium::Dp(unsigned int pair_index) {
        check(pair_index);
        return _Dp[pair_index];
    }

    double LinkageDisequilibrium::r(unsigned int pair_index) {
        check(pair_index);
        return _r[pair_index];
    }

    double LinkageDisequilibrium::r2(unsigned int pair_index) {
        check(pair_index);
        return _r2[pair_index];
    }

    unsigned int LinkageDisequilibrium::site1(unsigned int pair_index) {
        check(pair_index);
        return _site1[pair_index];
    }

    unsigned int LinkageDisequilibrium::site2(unsigned int pair_index) {
        check(pair_index);
        return _site2[pair_index];
    }

    double LinkageDisequilibrium::correl() const {
        if (_n==0) return 0.;
        double r, a;
        _correl(_n, _d, _r2, r, a);
        return r;
    }


    //  important helper

    void LinkageDisequilibrium::_correl(unsigned int n, const int* x, const double* y, double& r, double& a) {
        double X=0.;    // sum, then average
        double Y=0.;    // sun, then average
        double SSDx=0.; // sum of squared deviations
        double SSDy=0.; // sum of squared deviations
        double SJD=0.;  // sum of joint deviations

        for (unsigned int i=0; i<n; i++) {
            X+=x[i];
            Y+=y[i];
        }
        X/=n;
        Y/=n;

        for (unsigned int i=0; i<n; i++) {
            SSDx+= (x[i]-X)*(x[i]-X);
            SSDy+= (y[i]-Y)*(y[i]-Y);
            SJD+= (y[i]-Y)*(x[i]-X);
        }

        // coefficients
        a = SJD/SSDx;                 // regression coefficient (regression of y by x)
        r = SJD/(sqrt(SSDx)*sqrt(SSDy));    // correlation coefficient
    }
    
    
    // Hudson and Kaplan's minimal number of recombination events
    
    unsigned int LinkageDisequilibrium::Rmin(CharMatrix& data) const {

        // intervals between diallelic sites
        
        std::vector<Interval> intervals;
        unsigned int last=0;
        for (unsigned int i=0; i<v_S; i++) {
            if (v_sites[i]->numberOfAlleles()!=2) continue;
            for (unsigned int j=i+1; j<v_S; j++) {
                if (v_sites[j]->numberOfAlleles()!=2) continue;
                Interval interval(i,j);
                intervals.push_back(interval);
                last = j;
            }
        }

        // erases all intervals that have four alleles
        
        for (std::vector<Interval>::iterator it=intervals.begin(); 
                                            it!=intervals.end(); ++it) {
                                                
            std::vector<unsigned int> alleles(4);
            unsigned int ca = 0;
                    
                    // we try to tolerate missing data
                    
            for (unsigned int i=0; i<data.numberOfSequences(); i++) {

                if (data.populationLabel(i)==999) continue;
                
                char c1 = toupper(data.character(i, it->a()));
                char c2 = toupper(data.character(i, it->b()));

                if ((c1!='A' && c1!='C' && c1!='G' && c1!='T') ||
                    (c2!='A' && c2!='C' && c2!='G' && c2!='T')) continue;

                // try to find the game in previous (ca) gametes
                
                bool found=false;
                for (unsigned int j=0; j<ca; j++) {
                    if (c1 == toupper(data.character(alleles[j], it->a()))
                     && c2 == toupper(data.character(alleles[j], it->b()))) {
                         found=true;
                         break;
                     }
                }
                
                if (!found) {
                    alleles[ca] = i;
                    ca++;
                }
                
                if (ca==4) {
                    break;
                }
            }

            // if not 4 gametes, removes the interval

            if (ca<4) it->set_false();
        }

        // removes nested intervals (removes outer one)
        
        for(unsigned int i=0; i<intervals.size(); i++) {

            if (!intervals[i].good()) continue;

            for (unsigned int j=i+1; j<intervals.size(); j++) {
            
                if (!intervals[j].good()) continue;

                if (intervals[i].b() < intervals[j].a()) {
                    break;
                }

                if (intervals[i].b() >= intervals[j].b()) {
                    intervals[i].set_false();
                    break;
                }
            }
        }

        // adds a dummy interval at the end of the list 
        
        Interval dummy(last,last);
        intervals.push_back(dummy);
        
        // removes overlapping intervals
        
        for (unsigned int i=0; i<intervals.size(); i++) {
            
            if (!intervals[i].good()) continue;

            for (unsigned int j=i+1; j<intervals.size(); j++) {
                
                if (!intervals[j].good()) continue;

                if (intervals[j].a() < intervals[i].b()) {
                    intervals[j].set_false();
                }

            }
            
        }

        unsigned int c=0;
        for (std::vector<Interval>::iterator it=intervals.begin(); 
                                            it!=intervals.end(); ++it) {
            if (it->good()) c++;
        }

        return c - 1;
    }
    
    
    // RMIN HELPER

    LinkageDisequilibrium::Interval::Interval(unsigned int A, unsigned int B) {
        _a = A;
        _b = B;
        _good = true;
    }
    
    bool LinkageDisequilibrium::Interval::good() const {
        return _good;
    }

    void LinkageDisequilibrium::Interval::set_false() {
        _good = false;
    }
    
    unsigned int LinkageDisequilibrium::Interval::a() const {
        return _a;
    }
    
    unsigned int LinkageDisequilibrium::Interval::b() const {
        return _b;
    }


}
