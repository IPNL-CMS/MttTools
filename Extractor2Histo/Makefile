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
ROOTLIBS   = $(shell root-config --libs) -lMinuit -lRooFitCore -lRooFit

CXX           = g++
CXXFLAGS	    = -O2 -Wall -Werror -fPIC
#CXXFLAGS	    = -g -Wall -Werror -fPIC
LD			      = g++
LDFLAGS		    = -O2 $(shell root-config --ldflags) -Wl,--as-needed
SOFLAGS		    = 
AR            = ar
ARFLAGS       = -cq

LDDIR         = $(shell echo "-L" $$LD_LIBRARY_PATH | sed -e 's/:/ -L/g') -L$(ROOFIT_LIBDIR)/ -L../lib/

CXXFLAGS	   += $(ROOTCFLAGS) $(INCLUDES) -I.
LIBS  		    = $(LDDIR) $(ROOTLIBS)# -lAnalysis -lTreeAnalyzer_PhotonJet
GLIBS	    	  = $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard *.$(SrcSuf))
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

###

#all: drawPhotonJet_2bkg drawPhotonJetExtrap draw_all_methods_vs_pt draw_all_methods_vs_npv draw_crosschecks
all: process_all process_all_old

clean:
	@rm *.o;

process_all: Extractor2Histos.o process_all.o
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

process_all_old: Extractor2Histos_old.o process_all_old.o
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

process_all_old_semie: Extractor2Histos_old.o process_all_old_semie.o
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

process_all_old_semimu: Extractor2Histos_old.o process_all_old_semimu.o
	$(LD) $(LIBS) $(SOFLAGS) $(LDFLAGS) $+ -o $@

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
