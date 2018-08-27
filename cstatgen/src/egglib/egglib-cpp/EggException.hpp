/*
    Copyright 2009 Stéphane De Mita, Mathieu Siol

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

#ifndef EGGLIB_EGGEXCEPTION_HPP
#define EGGLIB_EGGEXCEPTION_HPP

#include <string>
#include <exception>

namespace egglib {

    /** \brief Base exception type for errors occurring in this library
     *
     * \ingroup core
     *
     */
    class EggException : public std::exception {
        public:
            /// Constructor with empty error message
            EggException();
            /// Creates the exception
            EggException(const char* message);
            /// Destructor
            ~EggException() throw() {}
            /// Gets error message
            virtual const char* what() const throw();
            
      protected:
            std::string message;

    };

  
    /** \brief Exception type for memory errors
     *
     * \ingroup core
     *
     */
    class EggMemoryError : public EggException {
        public:
            /// Creates the exception
            EggMemoryError();
            /// Destructor
            ~EggMemoryError() throw() {}
    };


    /** \brief Exception type for argument value errors
     *
     * \ingroup core
     *
     */
    class EggArgumentValueError : public EggException {
        public:
            /// Creates the exception
            EggArgumentValueError(const char* m );
            /// Destructor
            ~EggArgumentValueError() throw() {}
    };
    

    /** \brief Exception type for runtime errors
     * 
     * Runtime error definition is rather large. Includes bugs as well
     * as logical errors.
     *
     * \ingroup core
     *
     */
    class EggRuntimeError : public EggException {
        public:
            /// Creates the exception
            EggRuntimeError(const char* m );
            /// Destructor
            ~EggRuntimeError() throw() {}
    };


    /** \brief Exception type for file/string formatting errors
     *
     * \ingroup core
     *
     */
    class EggFormatError : public EggException {
        public:
            /// Creates the exception
            EggFormatError(const char* fileName, const char* expectedFormat, const char* m);
            /// Destructor
            ~EggFormatError() throw() {}
            /// Gets the file name
            std::string fileName() const;
            /// Gets the expected format
            std::string expectedFormat() const;
            /// Formats a longer string
            virtual const char* what_more() const;
            
        protected:
            std::string fname;
            std::string eformat;
    };


    /** \brief Exception type for errors while opening a file
     * 
     * \ingroup core
     *
     */
    class EggOpenFileError : public EggException {
        public:
            /// Creates the exception
            EggOpenFileError(const char* fileName );
            /// Destructor
            ~EggOpenFileError() throw() {}
    };
    
    
    /** \brief Exception type for unaligned sequences
     * 
     * \ingroup core
     * 
     */
    class EggUnalignedError : public EggException {
        public:
           /** \brief Creates the exception
            * 
            */
            EggUnalignedError();
            
           /** \brief Destructor
            * 
            */
            ~EggUnalignedError() throw() {}
    };

   /** \brief Exception type for invalid character
    * 
    * \ingroup core
    * 
    */
    class EggInvalidCharacterError : public EggException {
        public:
           /** \brief Creates the exception
            * 
            */
            EggInvalidCharacterError(char c, unsigned int seqIndex, unsigned int posIndex);
            
           /** \brief Destructor
            * 
            */
            ~EggInvalidCharacterError() throw() {}
    };

}



#endif
