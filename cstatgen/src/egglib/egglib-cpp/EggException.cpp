/*
    Copyright 2009,2011 Stéphane De Mita, Mathieu Siol

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
#include <sstream>

namespace egglib {

    EggException::EggException() {
        message= "";
    }

    EggException::EggException(const char* m) {
        message = "EggLib exception: ";
        message+= m;
    }
            
    const char* EggException::what() const  throw() {
        return message.c_str();
    }

    EggMemoryError::EggMemoryError() {
        message = "EggLib exception: error while allocating memory";
    }


    EggArgumentValueError::EggArgumentValueError(const char* m) {
        message = "EggLib exception: invalid value for function argument: ";
        message+= m;
    }

    EggRuntimeError::EggRuntimeError(const char* m) {
        message = "EggLib exception: runtime error: ";
        message +=  m;
    }

    EggFormatError::EggFormatError(const char* fileName, const char* expectedFormat, const char* m) {
        message = m;
        eformat = expectedFormat;
        fname = fileName;
    }

    std::string EggFormatError::fileName() const {
        return fname;
    }
    
    std::string EggFormatError::expectedFormat() const {
        return eformat;
    }
    
    const char* EggFormatError::what_more() const {
        std::string flap = "EggLib Exception: formatting error: data from ";
        flap += fname;
        flap += " is no valid ";
        flap += eformat;
        flap += ": ";
        flap += message;
        return flap.c_str();
    }

    EggOpenFileError::EggOpenFileError(const char* fileName ) {
        message = "EggLib exception: error while opening this file: ";
        message+= fileName;
    }

    EggUnalignedError::EggUnalignedError() {
        message = "Sequence doesn't match the alignment length";
    }

    EggInvalidCharacterError::EggInvalidCharacterError(char c, unsigned int seqIndex, unsigned int posIndex) {
        std::ostringstream stream;
        stream << "EggLib exception: invalid character found - ";
        stream << c;
        stream << " found at position " << posIndex << " of sequence number " << seqIndex;
        message = stream.str();
    }

}

