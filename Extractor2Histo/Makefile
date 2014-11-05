ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o
HeadSuf       = h

#ROOFIT_INCLUDE := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep INCLUDE= | sed 's|INCLUDE=||')

#ROOFIT_LIBDIR := $(shell cd $(CMSSW_BASE); scram tool info roofitcore | grep LIBDIR= | sed 's|LIBDIR=||')

#INCLUDES = -I$(ROOFIT_INCLUDE)/ -I$(CMSSW_BASE)/src/JetMETCorrections/GammaJetFilter/bin/ -I$(boost_header_LOC_INCLUDE)
#INCLUDES = -I$(shell root-config --incdir)

ROOTSYS  ?= ERROR_RootSysIsNotDefined

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --noldflags --libs)

CXX           = g++
CXXFLAGS	    = -O2 -g -Wall -fPIC --std=c++0x
LD			      = g++
LDDIR         = -L$(shell root-config --libdir) -L../external/lib -L$(BOOST_ROOT)/lib/ -L../SelectionMVA/lib -L../BTag/lib
LDFLAGS		    = -fPIC $(shell root-config --ldflags) $(LDDIR)
SOFLAGS		    = 
AR            = ar
ARFLAGS       = -cq

CXXFLAGS	   += $(ROOTCFLAGS) $(INCLUDES) -I. -I.. -I../SelectionMVA/include -I../BTag/include -I../external/ -I../include/ -I../external/include/ -I$(shell echo $(BOOST_ROOT))/include
LIBS  		    = $(ROOTLIBS) -lboost_filesystem -lboost_regex -lGenVector -lTMVA
GLIBS	    	  = $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard *.$(SrcSuf))
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

###

all: extractorToHisto substractHistos createSystHistograms

clean:
	@rm *.o;
	@rm *.d;
	@rm ../PUReweighting/*.o;
	@rm ../external/*.o;

extractorToHisto: ../external/TopTriggerEfficiencyProvider.o Extractor2Histos.o GaussianProfile.o ../PUReweighting/PUReweighter.o ../SelectionMVA/lib/libMVAReader.a ../SelectionMVA/lib/libBDTCuts.a
	$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@ -Wl,-Bstatic -lMVAReader -lBDTCuts -lBTagUtils -lyaml-cpp -Wl,-Bdynamic $(LIBS)

substractHistos: substractHistos.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@ -Wl,-Bstatic -lyaml-cpp -Wl,-Bdynamic $(LIBS)

createSystHistograms: createSystHistograms.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@ -Wl,-Bstatic -lyaml-cpp -Wl,-Bdynamic $(LIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Make the dependencies
%.d: %.cpp
	@echo "Generating dependencies for $<"
	@set -e; $(CXX) -M $(CXXFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

ifneq ($(MAKECMDGOALS), clean) 
-include $(DEPENDS) 
endif
