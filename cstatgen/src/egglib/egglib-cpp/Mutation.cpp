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


#include "Mutation.hpp"
#include "EggException.hpp"
#include <cstdlib>


namespace egglib {

    Mutation::Mutation() {
        init();
    }


    void Mutation::init() {
        //age = 0.;
        position = 0.;
        segmentIndex = 0;
        actualSiteIndex = 0;
        //edge = NULL;
    }

}
