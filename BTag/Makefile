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
LDDIR         = -L$(shell root-config --libdir) -L../external/lib -L$(BOOST_ROOT)/lib/
LDFLAGS		    = -fPIC $(shell root-config --ldflags) $(LDDIR)
SOFLAGS		    = 
AR            = ar
ARFLAGS       = -cq

CXXFLAGS	   += $(ROOTCFLAGS) $(INCLUDES) -Iinclude -I../include/ -I../external/include/ -I$(shell echo $(BOOST_ROOT))/include
LIBS  		    = $(ROOTLIBS) -lboost_filesystem -lboost_regex
GLIBS	    	  = $(ROOTGLIBS)
#------------------------------------------------------------------------------
SOURCES		= $(wildcard src/*.$(SrcSuf))
OBJECTS		= $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
DEPENDS		= $(SOURCES:.$(SrcSuf)=.d)
SOBJECTS	= $(SOURCES:.$(SrcSuf)=.$(DllSuf))

.SUFFIXES: .$(SrcSuf) .$(ObjSuf)

BTAG_UTILS_OBJECTS = src/BTagUtils.o

###

all: createBTaggingEffHistos libBTagUtils.a

clean:
	@rm  -f *.o $(OBJECTS);
	@rm -f *.d $(DEPENDS);
	@rm -f lib/*
	@rm -f ../PUReweighting/*.o;

createBTaggingEffHistos: createBTaggingEffHistos.o
	$(LD) $(SOFLAGS) $(LDFLAGS) $+ -o $@ $(LIBS)

libBTagUtils.a: $(BTAG_UTILS_OBJECTS)
	@echo "Generating libBTagUtils.a"
	@mkdir -p lib
	@ar cru lib/libBTagUtils.a $(BTAG_UTILS_OBJECTS)
	@ranlib lib/libBTagUtils.a

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
