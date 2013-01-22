# Which libraries to include
USERLIBS= $(shell root-config --libs) -lRooFitCore -lRooFit -lFoam -lMinuit -lMathMore
#USERLIBS= $(shell root-config --libs) -lRooFit

ifdef CMSSW_RELEASE_BASE
# Where to search for the libraries
UNPATCHED_CMSSW_RELEASE_BASE=$(shell cat $(CMSSW_RELEASE_BASE)/config/toolbox/$(SCRAM_ARCH)/tools/selected/cmssw.xml | grep 'name=\"CMSSW_BASE\"' | awk -F' ' '{print $$3};' | awk -F'\"' '{print $$2;}')
USERLIBDIRS= -L$(CMSSW_BASE)/lib/$(SCRAM_ARCH) \
	-L$(CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) \
	-L$(UNPATCHED_CMSSW_RELEASE_BASE)/lib/$(SCRAM_ARCH) -L$(ROOFITSYS)/lib -L.

INCDIRS = -I$(shell root-config --incdir) -I$(ROOFITSYS)/include -I$(shell pwd)/include/
else
USERLIBDIRS = $(root-config --libdir) -L.

# Where to search for the include files.
INCDIRS = -I$(shell root-config --incdir) -I$(shell pwd)/include/
endif

CC = gcc
CXX = g++
LD = g++
CFLAGS = -fPIC $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags) -Wl,--as-needed
#LDFLAGS = $(shell root-config --ldflags)

# Here's where you'd want to define any #define with -DVariableName
# (which means you can now use '#ifdef VariableName' preprocessor
# lines in your source code).
CFLAGS += $(INCDIRS) -std=c++0x -g -Wall -Wshadow -Wpointer-arith -Wcast-qual -Wunused-variable -Wextra

# Everything it needs to know about libraries
DYNAMIC_LIBS = $(USERLIBDIRS) $(USERLIBS) -lboost_regex
LIBS =  -Wl,-Bstatic -lUtils -Wl,-Bdynamic $(DYNAMIC_LIBS)

# Define what to compile as well as the necessary dependencies and
# object names.
SRCS =  $(wildcard *.cc) 
DEPS =  $(patsubst %.cc, %.d, $(wildcard *.cc)) 
OBJS =  $(patsubst %.cc, %.o, $(wildcard *.cc))

all: fritSignal fitMtt computeEff computeCrystalBallSyst computeJECSystWithToys computeJECSyst computeBkgSyst computeSigmaRef runLikelihoodScan treatToyStuff extractLimits drawLimitCurve computeKeysPdfSyst interpolatePdf generateToyDataset

# Make the executable
fitMtt: fitMtt.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) fitMtt.o $(LIBS) $(EXPLLINKLIBS) -o fitMtt

generateToyDataset: generateToyDataset.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) generateToyDataset.o $(LIBS) $(EXPLLINKLIBS) -o generateToyDataset

interpolatePdf: interpolatePdf.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) interpolatePdf.o $(LIBS) $(EXPLLINKLIBS) -o interpolatePdf

fritSignal: fritSignal.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) fritSignal.o $(LIBS) $(EXPLLINKLIBS) -o fritSignal

computeEff: computeEff.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) computeEff.o $(LIBS) $(EXPLLINKLIBS) -o computeEff

computeCrystalBallSyst: computeCrystalBallSyst.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) computeCrystalBallSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeCrystalBallSyst

computeKeysPdfSyst: computeKeysPdfSyst.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) computeKeysPdfSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeKeysPdfSyst

computeJECSyst: computeJECSyst.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) computeJECSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeJECSyst

computeJECSystWithToys: computeJECSystWithToys.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) computeJECSystWithToys.o $(LIBS) $(EXPLLINKLIBS) -o computeJECSystWithToys

computeBkgSyst: computeBkgSyst.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) computeBkgSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeBkgSyst

computeSigmaRef: computeSigmaRef.o libUtils
	$(LD) $(LDFLAGS) computeSigmaRef.o $(LIBS) $(EXPLLINKLIBS) -o computeSigmaRef

runLikelihoodScan: runLikelihoodScan.o lib/jsoncpp.o
	$(LD) $(LDFLAGS) $^ $(DYNAMIC_LIBS) $(EXPLLINKLIBS) -o runLikelihoodScan

treatToyStuff: treatToyStuff.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) treatToyStuff.o $(LIBS) $(EXPLLINKLIBS) -o treatToyStuff

extractLimits: extractLimits.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) extractLimits.o $(LIBS) $(EXPLLINKLIBS) -o extractLimits

drawLimitCurve: drawLimitCurve.o lib/jsoncpp.o libUtils
	$(LD) $(LDFLAGS) drawLimitCurve.o $(LIBS) $(EXPLLINKLIBS) -o drawLimitCurve

libUtils: Utils.o lib/jsoncpp.o
	$(LD) -shared -fvisibility=default $(LDFLAGS) $^ $(EXPLLINKLIBS) -o libUtils.so
	-@rm libUtils.a 2>/dev/null || true
	ar -cq libUtils.a $^

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
	@rm *.so *.a

# If we aren't doing a 'gmake clean' then make the dependencies
ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

