/*
    Copyright 2009-2010 Stéphane De Mita, Mathieu Siol

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

#include "EggException.hpp"
#include "Random.hpp"
#include "Controller.hpp"
#include "Population.hpp"
#include <cmath>


namespace egglib {

    /**********************************************************************/

    Controller::Controller() {
        random=NULL;
        current.reset(&paramSet);
        arg.set(&current, paramSet.numberOfSegments());
    }

                
    /**********************************************************************/

    Controller::~Controller() {
    }


    /**********************************************************************/

    Controller::Controller(const ParamSet* paramSet, Random* random) {
        this->initialParamSet = paramSet;
        this->paramSet = *paramSet;
        this->random = random;
        current.reset(&this->paramSet);
        arg.set(&current, this->paramSet.numberOfSegments());
        diploids();
    }
    

    /**********************************************************************/
    
    void Controller::reset() {
        paramSet = *initialParamSet;
        current.reset(&paramSet);
        arg.reset(&current);
        diploids();
    }


    /**********************************************************************/

    void Controller::diploids() {
        double s = paramSet.selfingRate();
        for (unsigned int i=0; i<current.numberOfPopulations(); i++) {
            unsigned int c=0;
            for (unsigned int j=0; j<paramSet.doubles(i); j++) {
                if (random->uniform()<s/(2.-s)) {
                    arg.coalescence(0., i, c, c+1);
                }
                else c+=2;
            }
        }
    }


    /**********************************************************************/

    unsigned int Controller::step() {

        /******************************************************************\
        |******* Welcome to the heart of the coalescent simulator ! *******|
        \******************************************************************/

        // Safety escape - the function does nothing if there is nothing to coalesce
        if (current.totalNumberOfLineages()<2) return current.totalNumberOfLineages();

        // These two variables stores which event is supposed to occur next
        // and when it is planned to occur (as an offset from now, whenever
        // "now" is.
        // The code is
        //      0 no event (this will be treated as an error)
        //      C coalescence
        //      R recombination
        //      M migration
        //      D demographic change
        char next = '0';
        double time = -1;

        // Gets coalescent time
        double coalescenceTime = 0.;
        unsigned int coalescencePopulation = 0;
        getCoalescenceTime(coalescenceTime, coalescencePopulation);
        if (coalescenceTime!=-1) {
            next = 'C';
            time = coalescenceTime;
        }

        // Gets recombination time
        double recombinationTime = getRecombinationTime();
        if (recombinationTime!=-1 && (time==-1 || recombinationTime<time)) {
            next = 'R';
            time = recombinationTime;
        }

        // Gets migration time
        double migrationParameter=0.; // used later to draw the population (passed to migrate())
        double migrationTime = getMigrationTime(migrationParameter);
        if (migrationTime!=-1 && (time==-1 || migrationTime<time)) {
            next = 'M';
            time = migrationTime;
        }

        // Checks whether a demographic event occurs before the smallest
        // found time (the date must be converted to a time (from now)
        double changeTime = paramSet.nextChangeDate();
        if (changeTime!=-1) {
            changeTime-=arg.time();
            if (time==-1 || changeTime<time) {
                next = 'D';
                time = changeTime;
            }
        }


        // checks that we are not engaged in an infinite process - it
        // doesn't matter if a recombination time is planned
        if (coalescenceTime==-1 && migrationTime==-1 && changeTime==-1) {
            throw EggRuntimeError("cannot coalesce: the current parameters don't allow to terminate the simulation");
        }

        // performs the event
        switch (next) {
            
            // case of a coalescence
            case 'C':
                arg.coalescence(coalescenceTime, coalescencePopulation, random);
                break;
                
            // case of a migration
            case 'M':
                arg.addTime(migrationTime);
                migrate(migrationParameter);
                break;
                
            // case of a recombination
            case 'R':
                arg.recombination(recombinationTime, random);
                break;
                
            // case of a demographic change
            case 'D':
                arg.addTime(changeTime);
                paramSet.nextChangeDo(this);
                break;

            // safety check
            case '0':
            default:
                // this is not supposed to happen
                throw EggRuntimeError("a bug occurred within Controller::step()"); 
        }

        // terminates 
        
        return current.totalNumberOfLineages();
    }


    /**********************************************************************/

    Arg* Controller::getArg() {
        return &arg;
    }


    /**********************************************************************/

    void Controller::bottleneck(unsigned int populationIndex, double strength) {
        double localTime= 0.;
        double increment= -1;
        while (true) {
            increment = getCoalescenceTimeForPopulation(populationIndex);
            if (increment==-1) break; // no possible coalescence is interpreted as n=1 (therefore non-bug)
            if (localTime+increment>strength) break;
            arg.coalescence(0., populationIndex, random);
            localTime+=increment;
        }
    }


    /**********************************************************************/

    void Controller::moveAllLineages(unsigned int source, unsigned int dest) {
        while (current.population(source)->numberOfLineages() > 0) {
            Edge* edge = current.population(source)->extractByIndex(0);
            current.population(dest)->push(edge);
        }
    }


    /**********************************************************************/

    void Controller::moveSomeLineages(unsigned int source, unsigned int dest, double probability) {
        unsigned int i=0;
        while (i < current.population(source)->numberOfLineages()) {
            if (random->uniform()<probability) {
                Edge* edge = current.population(source)->extractByIndex(i);
                current.population(dest)->push(edge);
            }
            else i++;
        }
    }


    /**********************************************************************/

    void Controller::addPopulation() {
        current.addPopulation();
    }


    /**********************************************************************/

    double Controller::getMigrationTime(double& migrationParameterDestination) {
        // if not enough populations, bails out immediately
        if (current.numberOfPopulations()<2) return -1;
        
        // mig will be the sum of the diagonal weighted by the number of individuals in each pop
        migrationParameterDestination = 0.;
        for (unsigned int i=0; i<current.numberOfPopulations(); i++) {
            migrationParameterDestination+= current.population(i)->numberOfLineages() * paramSet.pairwiseMigrationRate(i,i);
        }

        
        // if not possibility for migration
        if (migrationParameterDestination<0.000000000001) return -1;
        
        // draw a migration time
        return random->erand(1./migrationParameterDestination);
    }


    /**********************************************************************/

    void Controller::getCoalescenceTime(double& destTime, unsigned int& destPopIndex) {
        // sets the default value to the return variables
        destTime = -1;
        destPopIndex = 0;
        
        for (unsigned int i=0; i<current.numberOfPopulations(); i++) {

            // draw a coalescent time for pop i
            double t = getCoalescenceTimeForPopulation(i);

            // if this time is shorter, set it
            if (t!=-1 && ( destTime==-1 || t<destTime)) {
                destTime = t;
                destPopIndex = i;
            }
        }
    }


    /**********************************************************************/

    double Controller::getCoalescenceTimeForPopulation(unsigned int populationIndex) {
        // computes the basic coalescence time (if enough lineages)
        unsigned int n = current.population(populationIndex)->numberOfLineages();
        if (n<2) return -1;
        double s = paramSet.selfingRate();
        double N = paramSet.populationSize(populationIndex);
        if (N==0.) return -1;

        double Tcoal = -1;
        double expectation = N*(2.-s) / (2.*n*(n-1));
        
        double A = paramSet.growthRate(populationIndex);
        if (A==0.) {
            Tcoal = random->erand(expectation);
        }
                
        // applies a correction if the population is exponentially growing/shrinking
        else {
            double T0 = paramSet.dateOfLastChange(populationIndex);
            double T = arg.time();
		     double arg  = 1.+A*exp(-A*(T-T0))*random->erand(expectation);
		     if (arg > 0) Tcoal = log(arg)/A;
             else Tcoal = -1;

        }

        return Tcoal;
    }


    /**********************************************************************/

    double Controller::getRecombinationTime() const {
        if (paramSet.recombinationRate()==0.) return -1;
        
        double R = paramSet.recombinationRate();
        double s = paramSet.selfingRate();
        R *= 1-s/(2.-s);   // correction for probability that the recombined
                           // fragments coalesce back immediately
        if (R==0.) return -1;


        R *= 1. * current.efficientNumberOfLineages() / (paramSet.numberOfSegments()-1);
        R=1./R;
        return random->erand(R);
    }


    /**********************************************************************/

    void Controller::migrate(double mig) {

        // picks a pair of population by drawing a position in the sum
        double X = random->uniform() * mig;
        double acc = 0.;
        for (unsigned int pop1=0; pop1<paramSet.numberOfPopulations(); pop1++) {
            for (unsigned int pop2=0; pop2<paramSet.numberOfPopulations(); pop2++) {
                if (pop1==pop2) continue;
                acc+= current.population(pop1)->numberOfLineages() *  paramSet.pairwiseMigrationRate(pop1, pop2);
                // operates the migration (the lineage is picked by Population)
                if (X<acc) {
                    Edge* edge = current.population(pop1)->extractRandomly(random);
                    current.population(pop2)->push(edge);
                    return;
                }
            
            }
        }
         
        throw EggRuntimeError("hole in Controller::migrate()");
    }

}
