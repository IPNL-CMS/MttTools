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

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXX           = g++
CXXFLAGS	    = -O2 -g -Wall -fPIC -std=c++0x
LD			      = g++
LDFLAGS		    = $(shell root-config --ldflags)
SOFLAGS		    = 
AR            = ar
ARFLAGS       = -cq

LDDIR         = $(shell echo "-L" $$LD_LIBRARY_PATH | sed -e 's/:/ -L/g') -L$(ROOFIT_LIBDIR)/ -L../lib/

CXXFLAGS	   += $(ROOTCFLAGS) $(INCLUDES) -I. -I../include/ -I../external/ -I../ExtractorPostprocessing
LIBS  		    = $(LDDIR) $(ROOTLIBS) -lGenVector -lTMVA
GLIBS	    	  = $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard *.$(SrcSuf))
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

###

all: createMVATree trainMVA

clean:
	@rm *.o;
	@rm ../PUReweighting/*.o;

createMVATree: createMVATree.o ../PUReweighting/PUReweighter.o ../external/TopTriggerEfficiencyProvider.o ../ExtractorPostprocessing/ExtractorPostprocessing.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@ $(LIBS)

trainMVA: trainMVA.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@ $(LIBS)

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
