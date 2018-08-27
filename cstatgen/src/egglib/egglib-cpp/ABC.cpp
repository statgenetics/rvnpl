/*
    Copyright 2011 Stéphane De Mita
    
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


#include "ABC.hpp"
#include <cstdlib>
#include "EggException.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cmath>
#include <vector>
#include <algorithm>
#include "gsl/gsl_multifit.h"
#include "gsl/gsl_matrix.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_errno.h"


namespace egglib {

	void egghandler(const char* reason, const char* file, int line, int gsl_errno) {
		std::ostringstream stream;
		stream << "[GSL error] " << gsl_strerror(gsl_errno) << " - " << reason;
		throw EggRuntimeError(stream.str().c_str());
	}

    ABC::ABC() {
        init();
    }


    ABC::~ABC() {
        clear();
    }

    void ABC::init() {
        _sd = NULL;
        _obs = NULL;
        _nstats = 0;
        _nsam = 0;
        _threshold = -1;
        gsl_set_error_handler(&egghandler);
    }
    
    
    void ABC::clear() {
        if (_sd) free(_sd);
        if (_obs) free(_obs);
        _fnames.clear();
        _nparams.clear();
    }

     
    void ABC::number_of_statistics(unsigned int ns) {
        clear();
        init();

        _nstats = ns;
        
        if (ns==0) throw EggArgumentValueError("ABC: at least one stat is needed");

        _sd = (double*) malloc(ns * sizeof(double));
        if (!_sd) throw EggMemoryError();

        _obs = (double*) malloc(ns * sizeof(double));
        if (!_obs) throw EggMemoryError();

        for (unsigned int i=0; i<ns; i++) {
            _sd[i] = 0.;
        }
        
    }


    void ABC::add_fname(const char* fname, unsigned int np) {
        if (np==0) throw EggArgumentValueError("ABC: at least one parameter is needed");
        std::string str = fname;
        _fnames.push_back(str);
        _nparams.push_back(np);

    }


    double ABC::sd(unsigned int index) const {
        if (index >= _nstats) {
            throw EggArgumentValueError("invalid statistic index in ABC");
        }
        return _sd[index];
    }


    void ABC::obs(unsigned int index, double value) {
        if (index >= _nstats) {
            throw EggArgumentValueError("invalid statistic index in ABC");
        }
        _obs[index] = value;
    }



    void ABC::get_threshold(double tolerance) {
        
        // initial checking

        if (_nstats==0) throw EggArgumentValueError("ABC: at least one stat is needed");
        if (_fnames.size()==0) throw EggArgumentValueError("ABC: at least one data file is needed");
        if (_nsam!=0) throw EggArgumentValueError("ABC: cannot call threshold() several times without clearing");
        if (_threshold!=-1) throw EggArgumentValueError("ABC: cannot call threshold() several times withoutout clearing");

        // common variables

        std::string line;
        std::string token;
        double value;
        
        // temporary tables
        
        double *sum = (double*) malloc(_nstats*sizeof(double));
        if (!sum) throw EggMemoryError();

        double *sum2 = (double*) malloc(_nstats*sizeof(double));
        if (!sum2) throw EggMemoryError();
        
        for (unsigned int i=0; i<_nstats; i++) {
            sum[i] = 0;
            sum2[i] = 0;
        }

        // gets data for standard deviations
        
        for (unsigned int i=0; i<_fnames.size(); i++) {

            std::ifstream fstream(_fnames[i].c_str());
            if (!fstream.good()) {
                free(sum);
                throw EggOpenFileError(_fnames[i].c_str());
            }

            while (!fstream.eof()) {
                
                if (!fstream.good()) {
                    free(sum);
                    free(sum2);
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "expected a new line");
                }
                
                getline(fstream, line);
                std::istringstream sstream(line);
                            
                if (line.size() == 0) break;
                
                for (unsigned int j=0; j<_nparams[i]; j++) {
                    if (!sstream.good()) {
                        free(sum);
                        free(sum2);
                        throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                }

                if (!sstream.good()) {
                    free(sum);
                    free(sum2);
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                }
                
                sstream >> token;
                if (token != "#") {
                    free(sum);
                    free(sum2);
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line (invalid number of parameters?)");
                }

                for (unsigned int j=0; j<_nstats; j++) {
                    if (!sstream.good()) {
                        free(sum);
                        free(sum2);
                        throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                    value = atof(token.c_str());
                    sum[j] += value;
                    sum2[j] += (value*value);
                }
                
                if (sstream.good()) throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "line longer than expected");

                _nsam++;
            }

            fstream.close();
        }

        if (_nsam==0) throw EggArgumentValueError("ABC: no sample data found");
        
        // computes standard deviations
        
        double m, m2;
        
        for (unsigned int i=0; i<_nstats; i++) {
            m = sum[i] / _nsam;
            m2 = sum2[i] / _nsam;
            _sd[i] = sqrt(m2 - m*m);
        }
        
        // releases memory
        
        free(sum);
        free(sum2);
        
        // computes Euclidean distances

        std::vector<double> euclid(_nsam, 0.);
        unsigned int cnt = 0;

        for (unsigned int i=0; i<_fnames.size(); i++) {

            std::ifstream fstream(_fnames[i].c_str());
            if (!fstream.good()) throw EggOpenFileError(_fnames[i].c_str());

            while (!fstream.eof()) {
               
                if (!fstream.good()) throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "expected a new line");
                
                getline(fstream, line);
                std::istringstream sstream(line);
            
                if (line.size() == 0) break;

                if (cnt>=_nsam) throw EggArgumentValueError("ABC: trying to compute more euclidean distances than the number of sample data used for computed statistics standard deviatons");
                
                for (unsigned int j=0; j<_nparams[i]; j++) {
                    if (!sstream.good()) {
                        free(sum);
                        free(sum2);
                        throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                }

                if (!sstream.good()) {
                    free(sum);
                    free(sum2);
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                }
                
                sstream >> token;
                if (token != "#") {
                    free(sum);
                    free(sum2);
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line (invalid number of parameters?)");
                }

                euclid[cnt] = 0.;

                for (unsigned int j=0; j<_nstats; j++) {
                    if (!sstream.good()) throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                    sstream >> token;
                    value = atof(token.c_str());
                    if (_sd[j]>0) {
                        value /= _sd[j];
                        euclid[cnt] += ( (value-_obs[j]/_sd[j]) * 
                                        (value-_obs[j]/_sd[j]) );
                    }
                }
                    
                if (sstream.good()) throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "line longer than expected");

                euclid[cnt] = sqrt(euclid[cnt]);
                cnt++;
            }
        
            fstream.close();
                
        }

        if (cnt!=_nsam) throw EggRuntimeError("ABC: inconsistent number of data between two steps to data access");
    
		unsigned int position = (unsigned int)ceil(_nsam * tolerance);
		
		if (position<2) {
			throw EggArgumentValueError("tolerance is not large enough to catch at least two points");
		}

        std::sort(euclid.begin(), euclid.end());
        _threshold = euclid[position];

    }


    unsigned int ABC::rejection(const char* outfname, bool exportlabels) {
        
        if (_threshold==-1) throw EggArgumentValueError("ABC: the rejection threshold must have been computed before rejection step");

        unsigned int n = 0;
        unsigned int accept = 0;
        double euclid;
        double value;
        double w;
        std::string line;
        std::string token;

        std::ofstream outfstream(outfname);
        if (!outfstream.good()) throw EggOpenFileError(outfname);

        for (unsigned int i=0; i<_fnames.size(); i++) {

            std::ifstream infstream(_fnames[i].c_str());
            if (!infstream.good()) throw EggOpenFileError(_fnames[i].c_str());

            while (!infstream.eof()) {
                
                if (!infstream.good()) throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "expected a new line");
                
                getline(infstream, line);
                std::istringstream sstream(line);
            
                if (line.size() == 0) break;

                for (unsigned int j=0; j<_nparams[i]; j++) {
                    if (!sstream.good()) {
                        throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                    }
                    sstream >> token;
                }

                if (!sstream.good()) {
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                }
                
                sstream >> token;
                if (token != "#") {
                    throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line (invalid number of parameters?)");
                }

                euclid = 0.;

                for (unsigned int j=0; j<_nstats; j++) {
                    if (!sstream.good()) throw EggFormatError(_fnames[i].c_str(), "ABC sample data", "invalid line");
                    sstream >> token;
                    value = atof(token.c_str());
                    if (_sd[j]>0) {
                        value /= _sd[j];
                        euclid += ( (value-_obs[j]/_sd[j]) * 
                                    (value-_obs[j]/_sd[j]) );
                    }
                }
                
                euclid = sqrt(euclid);
                n++;
                if (euclid < _threshold) {
                    accept++;
                    w = 1-((euclid*euclid)/(_threshold*_threshold));
                    if (exportlabels) outfstream << "[" << i+1 << "] ";
                    outfstream << line << " # " << w << std::endl;
                }
            }

            infstream.close();

        }

        if (n != _nsam) {
            throw EggRuntimeError("number of samples do not match between ABC::get_sd() and ABC::threshold() calls");
        }


        outfstream.close();
        
        return accept;        
    }



    unsigned int ABC::regression(const char* infname, const char* outfname,
                                 TransformMode mode, const char* header) {

        // first pass to get the number of point and checks file consistency
        
        std::ifstream fstream(infname);
        if (!fstream.good()) throw EggOpenFileError(infname);

        std::string line;
        std::string token;

        unsigned int nparams = 0;
        unsigned int npoints = 0;
        
        while (!fstream.eof()) {
            
            if (!fstream.good()) throw EggFormatError(infname, "ABC sample data", "expected a new line");
            
            getline(fstream, line);
            if (line.size()==0) break;
            if (line[0]=='[') throw EggFormatError(infname, "ABC sample data", "model labels are not supported for regression");
        
            std::istringstream sstream(line);
            unsigned int i=0;
        
            while (1) {
                if (!sstream.good()) throw EggFormatError(infname, "ABC sample data", "invalid line structure  (3)");
                sstream >> token;
                if (token == "#") break;
                i++;
            }
            
            if (i==0) throw EggFormatError(infname, "ABC sample data", "no params?");
            if (nparams!=0 && i!=nparams) throw EggFormatError(infname, "ABC sample data", "inconsistent number of parameters");
            if (nparams==0) nparams = i;
            npoints++;
        }

        if (npoints==0) throw EggRuntimeError("cannot perform ABC: no simulations");

        fstream.close();

        // allocates memory
        
		gsl_multifit_linear_workspace* workspace =
                    gsl_multifit_linear_alloc(npoints, _nstats+1);  // +1 for intercept
		
		gsl_matrix* mstats = gsl_matrix_alloc(npoints, _nstats+1); // +1 for intercept
		gsl_vector* vweights = gsl_vector_alloc(npoints);
        gsl_matrix* mparams = gsl_matrix_alloc(npoints, nparams);
        gsl_vector* vparams = gsl_vector_alloc(npoints);
		gsl_vector* vcoefs = gsl_vector_alloc(_nstats+1); // +1 for intercept
		gsl_matrix* vcov = gsl_matrix_alloc(_nstats+1, _nstats+1);  // +1 for intercept

		if (!mstats || !workspace || !vweights ||
            !vparams || !mparams || !vcoefs || !vcov) {
            if (mstats) gsl_matrix_free(mstats);
            if (vweights) gsl_vector_free(vweights);
            if (mparams) gsl_matrix_free(mparams);
            if (vparams) gsl_vector_free(vparams);
            if (vcoefs) gsl_vector_free(vcoefs);
            if (vcov) gsl_matrix_free(vcov);
            if (workspace) gsl_multifit_linear_free(workspace);
            throw EggMemoryError();
        }

        // loads data

        fstream.open(infname);
        if (!fstream.good()) throw EggOpenFileError(infname);

        unsigned int point=0;
        
        while (!fstream.eof()) {
            if (!fstream.good()) throw EggFormatError(infname, "ABC sample data", "expected a new line");
            
            getline(fstream, line);
            if (line.size() == 0) break;
            
            std::istringstream sstream(line);
        
            // gets params
        
            for (unsigned int param=0; param<nparams; param++) {
                if (!sstream.good()) throw EggFormatError(infname, "ABC sample data", "invalid line structure  (6)");
                sstream >> token;
                gsl_matrix_set(mparams, point, param, atof(token.c_str()));
            }

            sstream >> token;
            if (token != "#") throw EggFormatError(infname, "ABC sample data", "invalid line structure  (7)");

            // gets stats

            gsl_matrix_set(mstats, point, 0, 1.); // for intercept
			
            for (unsigned int stat=0; stat<_nstats; stat++) {
                if (!sstream.good()) throw EggFormatError(infname, "ABC sample data", "invalid line structure  (8)");
                sstream >> token;
                if (_sd[stat]==0) throw EggRuntimeError("at least one statistic not variable in the local region (or ABC class not used properly)");
                gsl_matrix_set(mstats, point, stat+1, atof(token.c_str())/_sd[stat]);
            }
            sstream >> token;
            if (token != "#") throw EggFormatError(infname, "ABC sample data", "invalid line structure  (9)");

            // gets weigth
            
            if (!sstream.good()) throw EggFormatError(infname, "ABC sample data", "invalid line structure  (10)");
            sstream >> token;
            gsl_vector_set(vweights, point, atof(token.c_str()));

            point++;
            
        }
        
        if (point!=npoints) throw EggRuntimeError("inconsistent ABC file parsing...");

        fstream.close();

		// fits all parameters separately

        double chisq;
		double mean;
        
		for (unsigned int param=0; param<nparams; param++) {

			// loads the explained variable

            for (unsigned int point=0; point<npoints; point++) {
				gsl_vector_set(vparams, point,
                                gsl_matrix_get(mparams, point, param));
            }

            // performs transformation
            
            double min=0.;
            double max=999999999.;
            double x=0.;
            
            if (mode==TAN) {
                min = gsl_vector_get(vparams, 0);
                max = gsl_vector_get(vparams, 0);
                for (unsigned int point=1; point<npoints; point++) {
                    x = gsl_vector_get(vparams, point);
                    if (x<min) min=x;
                    if (x>max) max=x;
                }
                min -= 0.00000001;
                max += 0.00000001;
                for (unsigned int point=0; point<npoints; point++) {
                    x = gsl_vector_get(vparams, point);
                    x = -log(1./(tan( ((x-min)/(max-min))*(M_PI/2.) )));
                    gsl_vector_set(vparams, point, x);
                }
            }

            if (mode==LOG) {
                for (unsigned int point=0; point<npoints; point++) {
                    x = gsl_vector_get(vparams, point);
					if (x<=0.) throw EggRuntimeError("ABC transformation LOG not available for parameter values <= 0");
					gsl_vector_set(vparams, point, log(x));
                }
            }
            
			// performs regression

			gsl_multifit_wlinear(mstats, vweights, vparams,
									vcoefs, vcov, &chisq, workspace);

			// computes predicted mean

			mean = 1. * gsl_vector_get(vcoefs, 0);   // intercept
			for (unsigned int stat=0; stat<_nstats; stat++) {
				mean += _obs[stat]/_sd[stat] * gsl_vector_get(vcoefs, stat+1);
			}

			// computes residuals and predicted params (overwrites params table)

			for (unsigned int point=0; point<npoints; point++) {
                double observed = gsl_vector_get(vparams, point);
				double predicted = 0.;
				for (unsigned int stat=0; stat<_nstats+1; stat++) {
					predicted += gsl_matrix_get(mstats, point, stat)
											 * gsl_vector_get(vcoefs, stat);
				}

                x = mean + (observed - predicted);
                
                // untransforms immediately
                
                if (mode==LOG) x = exp(x);
                if (mode==TAN) x = min + (2./M_PI)*(max-min)*atan(exp(x));

                // and sets final value

                gsl_matrix_set(mparams, point, param, x);
			}

		}
        
        // exports data

        std::ofstream out(outfname);
        if (!out.good()) throw EggOpenFileError(outfname);

        if (strlen(header)) out << header << std::endl;

        for (unsigned int point=0; point<npoints; point++) {
            out << gsl_matrix_get(mparams, point, 0);
            for (unsigned int param=1; param<nparams; param++) {
                out << " " << gsl_matrix_get(mparams, point, param);
            }
            out << std::endl;
        }


        out.close();

        // frees memory

		gsl_multifit_linear_free(workspace);
		gsl_matrix_free(mstats);
		gsl_vector_free(vweights);
		gsl_matrix_free(mparams);
		gsl_vector_free(vparams);
		gsl_vector_free(vcoefs);
		gsl_matrix_free(vcov);


        return npoints;
        
    }


    unsigned int ABC::number_of_samples() const {
        return _nsam;
    }

    double ABC::threshold() const {
        return _threshold;
    }

}
