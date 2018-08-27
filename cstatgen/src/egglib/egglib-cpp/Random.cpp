/*
    Copyright 2008,2009,2012 Stéphane De Mita, Mathieu Siol, adapted from
    MStrat, developed by Charles-Edouard Coste, Thomas M. Bataillon,
    Mathieu Cotisson, Guy Decoux, Chistophe Rozale, Daniel J. Schoen
    and Jacques L. David.
    
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


#include "Random.hpp"
#include <cmath>
#include <ctime>
#include <cstdlib>



namespace egglib {

    Random::Random() {

        unsigned int now = 1254313777; // time(0) when this code was written
        int  secondsSinceCoded = time(0)-now; // seconds elapsed
        unsigned int uSecondsSinceCoded = (secondsSinceCoded<0)?-secondsSinceCoded:secondsSinceCoded; // supports accidental time travel
        _seed1 = uSecondsSinceCoded / 3600; // seed1: hours since coded
        _seed2 = uSecondsSinceCoded - _seed1 * 3600; // seed2: seconds of the current hour
        
        b_ncached = false;
        v_ncached = 0.;
    }

    Random::Random(double s1, double s2) {
        _seed1 = s1;
        _seed2 = s2;

        b_ncached = false;
        v_ncached = 0.;
    }

    double Random::erand(double esperance) {
        double tp = 0.0;
        while (tp==0.0) tp = uniform();
        return ( -(esperance)*log(tp));
    }
          
    unsigned int Random::irand(unsigned int x) {
        return (unsigned int) (uniform()*x);
    }  

    unsigned int Random::prand(double mean) {
        unsigned int i=0;
        double cumul;
        cumul= (-1/mean)*log(uniform());
        while (cumul<1) {
            cumul += (-1/mean)*log(uniform());
            i++;
        }
        return i;
    }

    unsigned int Random::grand(double param) {
        if (param==1.) return 1;
        double X = 1.-uniform();
        return (unsigned int) ceil(log(X)/log(1.-param));
    }

    double Random::uniform() {
        double z, r;
        int i;
    
        r = _seed1 / 53668;
        i = (int) r; // integer part of r
        r = (double) i; // restores double precision
        _seed1 = 40014 * (_seed1 - r * 53668) - r * 12211; // redefines first seed
        if (_seed1 < 0)  _seed1 += 2147483563;
        r = _seed2 / 52774;
        i = (int) r;
        r = (double) i;
        _seed2 = 40692 * (_seed2 - r * 52774) - r * 3791;
        if (_seed2 < 0) _seed2 += 2147483399;
        z = _seed1 - _seed2;
        if (z < 1) z += 2147483562;

        return (z * 4.656613e-10);
    }
    
    double Random::nrand() {
        if (b_ncached) {
            b_ncached = false;
            return v_ncached;
        }
        
        // polar form of the Box-Muller transformation
        // I would like the Ziggurat algorithm
        
        // implementation taken as is from http://www.taygeta.com/random/gaussian.html Nov 10th 2010
        
        float x1, x2, w, y1, y2;
 
        do {
            x1 = 2.0 * uniform() - 1.0;
            x2 = 2.0 * uniform() - 1.0;
            w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         y1 = x1 * w;
         y2 = x2 * w;

        // caches one value, return the other
        b_ncached = true;
        v_ncached = y2;
        return y1;
    }

    double Random::seed1() const { return _seed1; }
    
    double Random::seed2() const { return _seed2; }
    
    void Random::seed1(double seed) {
        _seed1 = seed;
        b_ncached = false;
        v_ncached = 0.;
    }
    
    void Random::seed2(double seed) {
        _seed2 = seed;
        b_ncached = false;
        v_ncached = 0.;
    }

}
