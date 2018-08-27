#AUTOMAKE_OPTIONS=foreign

AM_CPPFLAGS = -Wall -fPIC -O4 $(OPENMP_CFLAGS)
AM_LDFLAGS = $(OPENMP_CFLAGS)
lib_LIBRARIES = libegglib-cpp.a

libegglib_cpp_a_SOURCES = Align.cpp Container.cpp Convert.cpp \
	DataMatrix.cpp Fasta.cpp Ms.cpp Random.cpp EggException.cpp \
	Staden.cpp BaseDiversity.cpp Consensus.cpp FStatistics.cpp \
	HFStatistics.cpp HaplotypeDiversity.cpp LinkageDisequilibrium.cpp \
	MicrosatelliteDiversity.cpp NucleotideDiversity.cpp \
	SitePolymorphism.cpp Arg.cpp ChangeTypes.cpp Controller.cpp \
	Current.cpp Edge.cpp Mutation.cpp Mutator.cpp ParamSet.cpp \
	Population.cpp

library_include_HEADERS = Align.hpp Container.hpp Convert.hpp \
	CharMatrix.hpp DataMatrix.hpp Fasta.hpp Ms.hpp Random.hpp \
	EggException.hpp Staden.hpp BaseDiversity.hpp Consensus.hpp \
	FStatistics.hpp HFStatistics.hpp HaplotypeDiversity.hpp \
	LinkageDisequilibrium.hpp MicrosatelliteDiversity.hpp \
	NucleotideDiversity.hpp SitePolymorphism.hpp Arg.hpp \
	ChangeTypes.hpp Controller.hpp Current.hpp Edge.hpp Mutation.hpp \
	Mutator.hpp ParamSet.hpp Population.hpp \
	config.h

library_includedir=$(includedir)/egglib-cpp

if BPP_AVAILABLE
libegglib_cpp_a_SOURCES += BppDiversity.cpp 
library_include_HEADERS += BppDiversity.hpp
endif

if GSL_AVAILABLE
libegglib_cpp_a_SOURCES += ABC.cpp 
library_include_HEADERS += ABC.hpp
endif


