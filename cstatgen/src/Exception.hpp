// $File: Exception.hpp $
// $LastChangedDate:  $
// $Rev:  $
// Copyright (c) 2014, Gao Wang <ewanggao@gmail.com>
// GNU General Public License (http://www.gnu.org/licenses/gpl.html)

#ifndef _EXCEPTION_HPP_
#define _EXCEPTION_HPP_
namespace cstatgen {
/// exception handler. Exceptions will be passed to Python.
class Exception
{
public:
	/// constructor
	/// \param msg error message
	Exception(const std::string & msg) : m_msg(msg)
	{
	}


	/// return error message
	const char * message()
	{
		return m_msg.c_str();
	}


	virtual ~Exception()
	{
	};

private:
	/// error message
	std::string m_msg;
};

/// exception, thrown if out of memory
class StopIteration : public Exception
{
public:
	StopIteration(const std::string msg) : Exception(msg)
	{
	};
};


/// exception, thrown if index out of range
class IndexError : public Exception
{
public:
	IndexError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if value of range etc
class ValueError : public Exception
{
public:
	ValueError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if system error occurs
class SystemError : public Exception
{
public:
	SystemError(const std::string msg) : Exception(msg)
	{
	};
};

/// exception, thrown if a runtime error occurs
class RuntimeError : public Exception
{
public:
	RuntimeError(const std::string msg) : Exception(msg)
	{
	};
};
}
#endif
