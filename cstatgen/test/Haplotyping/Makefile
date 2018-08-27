CXX = g++ 
CC = gcc
UMICHLIB=../../src/umich
CORELIB=../../src
CFLAGS=-O3 -I./ -I$(CORELIB) -I $(UMICHLIB)/general -I$(UMICHLIB)/klib -I$(UMICHLIB)/libsrc -I$(UMICHLIB)/merlin -I$(UMICHLIB)/pdf -I$(UMICHLIB)/clusters -D_FILE_OFFSET_BITS=64 -D__ZLIB_AVAILABLE__ -Wall
CXXFLAGS= $(CFLAGS) -std=c++11
BINDIR = ./ 
MERLIN = demo.exe
EXECUTABLES = $(MERLIN)
# MERLIN File Set
MERLINBASE = $(UMICHLIB)/merlin/AssociationAnalysis $(UMICHLIB)/merlin/FastAssociation \
 $(UMICHLIB)/merlin/AnalysisTask $(UMICHLIB)/merlin/Conquer \
 $(UMICHLIB)/merlin/ConquerHaplotyping $(UMICHLIB)/merlin/DiseaseModel \
 $(UMICHLIB)/merlin/ParametricLikelihood $(UMICHLIB)/merlin/GenotypeInference \
 $(UMICHLIB)/merlin/Houdini \
 $(UMICHLIB)/merlin/KongAndCox $(UMICHLIB)/merlin/Manners $(UMICHLIB)/merlin/MerlinBitSet \
 $(UMICHLIB)/merlin/MerlinCluster $(UMICHLIB)/merlin/MerlinCore \
 $(UMICHLIB)/merlin/MerlinError $(UMICHLIB)/merlin/MerlinFamily $(UMICHLIB)/merlin/MerlinIBD \
 $(UMICHLIB)/merlin/InformationContent $(UMICHLIB)/merlin/MerlinCache $(UMICHLIB)/merlin/MerlinModel \
 $(UMICHLIB)/merlin/MerlinKinship $(UMICHLIB)/merlin/MerlinKinship15 \
 $(UMICHLIB)/merlin/MerlinHaplotype $(UMICHLIB)/merlin/MerlinMatrix $(UMICHLIB)/merlin/MerlinParameters \
 $(UMICHLIB)/merlin/MerlinPDF $(UMICHLIB)/merlin/MerlinSimulator $(UMICHLIB)/merlin/MerlinSimwalk2 \
 $(UMICHLIB)/merlin/MerlinSort $(UMICHLIB)/merlin/NPL-ASP $(UMICHLIB)/merlin/NPL-QTL \
 $(UMICHLIB)/merlin/Magic $(UMICHLIB)/merlin/Mantra $(UMICHLIB)/merlin/Parametric $(UMICHLIB)/merlin/QtlModel \
 $(UMICHLIB)/merlin/Tree \
 $(UMICHLIB)/merlin/TreeBasics $(UMICHLIB)/merlin/TreeIndex $(UMICHLIB)/merlin/TreeManager \
 $(UMICHLIB)/merlin/TreeInfo $(UMICHLIB)/merlin/TreeFlips $(UMICHLIB)/merlin/VarianceComponents
MERLINHDR = $(MERLINBASE:=.h) $(UMICHLIB)/merlin/TreeNode.h $(CORELIB)/HaplotypingUtils.hpp
MERLINSRC = $(MERLINBASE:=.cpp) demo.cpp $(CORELIB)/HaplotypingUtils.cpp
MERLINOBJ = $(MERLINSRC:.cpp=.o)

# Files for dealing with clustered markers
CLUSTERS = $(UMICHLIB)/clusters/HaploFamily $(UMICHLIB)/clusters/HaploGraph \
 $(UMICHLIB)/clusters/HaploSet $(UMICHLIB)/clusters/HaploTree $(UMICHLIB)/clusters/Likelihood \
 $(UMICHLIB)/clusters/SparseLikelihood $(UMICHLIB)/clusters/Unknown
CLUSTERCPP = $(CLUSTERS:=.cpp)
CLUSTERHDR = $(CLUSTERS:=.h)
CLUSTEROBJ = $(CLUSTERS:=.o)
# Utility Library File Set
LIBFILE = $(UMICHLIB)/general/lib-goncalo.a
LIBMAIN = $(UMICHLIB)/general/BasicHash $(UMICHLIB)/general/Error $(UMICHLIB)/general/FortranFormat \
 $(UMICHLIB)/general/GenotypeLists $(UMICHLIB)/general/InputFile $(UMICHLIB)/general/IntArray $(UMICHLIB)/general/Hash \
 $(UMICHLIB)/libsrc/LongArray $(UMICHLIB)/libsrc/Kinship $(UMICHLIB)/libsrc/KinshipX  $(UMICHLIB)/general/MapFunction \
 $(UMICHLIB)/libsrc/MathCholesky $(UMICHLIB)/libsrc/MathDeriv $(UMICHLIB)/libsrc/MathFloatVector \
 $(UMICHLIB)/libsrc/MathGenMin $(UMICHLIB)/libsrc/MathGold $(UMICHLIB)/general/MathMatrix $(UMICHLIB)/libsrc/MathStats \
 $(UMICHLIB)/libsrc/MathNormal $(UMICHLIB)/libsrc/MathSVD $(UMICHLIB)/general/MathVector \
 $(UMICHLIB)/general/MemoryInfo $(UMICHLIB)/general/MiniDeflate \
 $(UMICHLIB)/general/Parameters $(UMICHLIB)/general/Pedigree $(UMICHLIB)/general/PedigreeAlleleFreq \
 $(UMICHLIB)/general/PedigreeDescription $(UMICHLIB)/general/PedigreeFamily $(UMICHLIB)/general/PedigreeGlobals \
 $(UMICHLIB)/general/PedigreePerson $(UMICHLIB)/general/QuickIndex $(UMICHLIB)/general/Random $(UMICHLIB)/general/Sort \
 $(UMICHLIB)/general/StringArray $(UMICHLIB)/general/StringBasics $(UMICHLIB)/general/StringMap \
 $(UMICHLIB)/general/StringHash $(UMICHLIB)/libsrc/TraitTransformations \
 $(UMICHLIB)/klib/bgzf $(UMICHLIB)/klib/knetfile $(UMICHLIB)/general/PhoneHome \
 $(UMICHLIB)/general/GzipHeader $(UMICHLIB)/general/BgzfFileType $(UMICHLIB)/general/BgzfFileTypeRecovery \
 $(UMICHLIB)/general/GzipFileType $(UMICHLIB)/general/UncompressedFileType $(UMICHLIB)/general/FileType

LIBPED = $(UMICHLIB)/general/PedigreeLoader $(UMICHLIB)/general/PedigreeTwin $(UMICHLIB)/general/PedigreeTrim
LIBSRC = $(LIBMAIN:=.cpp) $(LIBPED:=.cpp)
LIBHDR = $(LIBMAIN:=.h) $(UMICHLIB)/general/Constant.h \
 $(UMICHLIB)/general/MathConstant.h $(UMICHLIB)/general/PedigreeAlleles.h $(UMICHLIB)/general/LongInt.h
LIBOBJ = $(LIBSRC:.cpp=.o)
# PDF Library File Sets
PDFLIB = $(UMICHLIB)/pdf/libpdf.a
PDFFILES = $(UMICHLIB)/pdf/PDF $(UMICHLIB)/pdf/PDFfont $(UMICHLIB)/pdf/PDFinfo $(UMICHLIB)/pdf/PDFpage \
 $(UMICHLIB)/pdf/PDFchartbasics $(UMICHLIB)/pdf/PDFchartbar $(UMICHLIB)/pdf/PDFlinechart \
 $(UMICHLIB)/pdf/PDFhistogram \
 $(UMICHLIB)/pdf/PDFchartaxis $(UMICHLIB)/pdf/PDFchartlegend $(UMICHLIB)/pdf/PDFchartmarker \
 $(UMICHLIB)/pdf/PDFchartline $(UMICHLIB)/pdf/PDFchartobject
PDFSRC = $(PDFFILES:=.cpp)
PDFHDR = $(PDFFILES:=.h)
PDFOBJ = $(PDFFILES:=.o)

# make everything
all : $(EXECUTABLES)

$(EXECUTABLES) : $(BINDIR)

$(BINDIR) :
	mkdir -p $(BINDIR)

# dependencies for executables
$(MERLIN) : $(LIBFILE) $(PDFLIB) $(MERLINOBJ) $(CLUSTEROBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(MERLINOBJ) $(CLUSTEROBJ) $(PDFLIB) $(LIBFILE) -lm -lz

$(LIBFILE) : $(LIBOBJ) $(LIBHDR)
	ar -cr $@ $(LIBOBJ)
	ranlib $@

$(PDFLIB) : $(PDFOBJ)
	ar -cr $@ $(PDFOBJ)
	ranlib $@

$(MERLINOBJ) : $(MERLINHDR) $(CLUSTERHDR) $(LIBHDR)

$(CLUSTEROBJ) : $(CLUSTERHDR) $(MERLINHDR) $(LIBHDR)

$(LIBOBJ) : $(LIBHDR)

$(PDFOBJ) : $(PDFHDR)

clean :
	rm -f ./*.o $(UMICHLIB)/*/*.a $(UMICHLIB)/*/*.o $(EXECUTABLES) $(CORELIB)/*.o

.c.o :
	$(CC) $(CFLAGS) -o $@ -c $*.c

#.cpp.X.o : 
#	$(CXX) $(CXXFLAGS) -o $@ -c $*.cpp -D__CHROMOSOME_X__

.cpp.o : 
	$(CXX) $(CXXFLAGS) -o $@ -c $*.cpp 

.SUFFIXES : .cpp .c .o .X.o $(SUFFIXES)
