# Which libraries to include
USERLIBS= $(shell root-config --libs) -lRooFitCore -lRooFit -lFoam -lMinuit -lMathMore
#USERLIBS= $(shell root-config --libs) -lRooFit

ifdef CMSSW_RELEASE_BASE
# Where to search for the libraries
UNPATCHED_CMSSW_RELEASE_BASE=$(shell cat $(CMSSW_RELEASE_BASE)/config/toolbox/$(SCRAM_ARCH)/tools/selected/cmssw.xml | grep 'name=\"CMSSW_BASE\"' | awk -F' ' '{print $$3};' | awk -F'\"' '{print $$2;}')
USERLIBDIRS= -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH) \
	-L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) \
	-L$(UNPATCHED_CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L$(ROOFITSYS)/lib

INCDIRS = -I$(shell root-config --incdir) -I$(ROOFITSYS)/include -I$(shell pwd)/include/
else
USERLIBDIRS = $(root-config --libdir)

# Where to search for the include files.
INCDIRS = -I$(shell root-config --incdir) -I$(shell pwd)/include/
endif

CC = gcc
CXX = g++
LD = g++
CFLAGS = $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags) -Wl,--as-needed
#LDFLAGS = $(shell root-config --ldflags)

# Here's where you'd want to define any #define with -DVariableName
# (which means you can now use '#ifdef VariableName' preprocessor
# lines in your source code).
CFLAGS += $(INCDIRS) -std=c++0x -g -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wunused-variable

# Everything it needs to know about libraries
LIBS =  $(USERLIBDIRS) $(USERLIBS) -lboost_regex

# Define what to compile as well as the necessary dependencies and
# object names.
SRCS =  $(wildcard *.cc) 
DEPS =  $(patsubst %.cc, %.d, $(wildcard *.cc)) 
OBJS =  $(patsubst %.cc, %.o, $(wildcard *.cc))

all: reduceDataset reduceDatasetFromData fritSignal fitMtt computeEff computeSignalSyst computeJECSyst computeBkgSyst computeSystPb computeSigmaRef runLikelihoodScan treatToyStuff extractLimits drawLimitCurve

# Make the executable
fitMtt: fitMtt.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o fitMtt

reduceDataset: reduceDataset.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o reduceDataset

reduceDatasetFromData: reduceDatasetFromData.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o reduceDatasetFromData

fritSignal: fritSignal.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o fritSignal

computeEff: computeEff.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o computeEff

computeSignalSyst: computeSignalSyst.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o computeSignalSyst

computeJECSyst: computeJECSyst.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o computeJECSyst

computeBkgSyst: computeBkgSyst.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o computeBkgSyst

computeSystPb: computeSystPb.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o computeSystPb

computeSigmaRef: computeSigmaRef.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o computeSigmaRef

runLikelihoodScan: runLikelihoodScan.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o runLikelihoodScan

treatToyStuff: treatToyStuff.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o treatToyStuff

extractLimits: extractLimits.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o extractLimits

drawLimitCurve: drawLimitCurve.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(LIBS) $(EXPLLINKLIBS) -o drawLimitCurve

# Make the objects
%.o: %.cc
	$(CXX) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CFLAGS) -c -o $@ $<

# Make the dependencies
%.d: %.cc
	@echo "Generating dependencies for $<"
	@set -e; $(CXX) -M $(CFLAGS) $< \
	| sed 's%\($*\)\.o[ :]*%\1.o $@ : %g' > $@; \
	[ -s $@ ] || rm -f $@

# For debugging purposes
echo:
	@echo For Debugging:
	@echo .
	@echo SRCS: $(SRCS)
	@echo .
	@echo INCLUDES: $(INCDIRS)
	@echo .
	@echo OBJECTS: $(OBJS)
	@echo .
	@echo DEPS: $(DEPS)
	@echo .
	@echo CFLAGS: $(CFLAGS)

# Clean everythingg
clean:
	@rm -f core* ${DEPS} ${OBJS} lib/*.o lib/*.d

# If we aren't doing a 'gmake clean' then make the dependencies
ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

