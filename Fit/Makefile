# Which libraries to include
USERLIBS= $(shell root-config --libs) -lRooFitCore -lRooFit -lFoam -lMinuit
#USERLIBS= $(shell root-config --libs) -lRooFit

USERLIBDIRS = $(root-config --libdir) -L$(BOOST_ROOT)/lib -L.

# Where to search for the include files.
INCDIRS = -I$(shell root-config --incdir) -I../external/include/ -I$(BOOST_ROOT)/include -I.

CC = gcc
CXX = g++
LD = g++
CFLAGS = -fPIC $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags) -L../external/lib -Wl,--as-needed
#LDFLAGS = $(shell root-config --ldflags)

# Here's where you'd want to define any #define with -DVariableName
# (which means you can now use '#ifdef VariableName' preprocessor
# lines in your source code).
CFLAGS += $(INCDIRS) -std=c++0x -O2 -g -Wall -Wshadow -Wpointer-arith -Wunused-variable -Wextra

# Everything it needs to know about libraries
DYNAMIC_LIBS = $(USERLIBDIRS) $(USERLIBS) -lboost_regex
LIBS =  -Wl,-Bstatic -ljsoncpp -lUtils -Wl,-Bdynamic $(DYNAMIC_LIBS)

# Define what to compile as well as the necessary dependencies and
# object names.
SRCS =  $(wildcard *.cc) 
DEPS =  $(patsubst %.cc, %.d, $(wildcard *.cc)) 
OBJS =  $(patsubst %.cc, %.o, $(wildcard *.cc))

all: fritSignal fitMtt createCombineWorkspace computeEff computeCrystalBallSyst computeJECSystWithToys computeJECSyst computeBkgSyst computeSigmaRef runLikelihoodScan treatToyStuff extractLimits drawLimitCurve computeKeysPdfSyst interpolatePdf generateToyDataset generateKeysPdfWorkspaces higgsSignalToPDF

# Make the executable
fitMtt: fitMtt.o libUtils.a
	$(LD) $(LDFLAGS) fitMtt.o $(LIBS) $(EXPLLINKLIBS) -o fitMtt

createCombineWorkspace: createCombineWorkspace.o libUtils.a
	$(LD) $(LDFLAGS) createCombineWorkspace.o $(LIBS) $(EXPLLINKLIBS) -o createCombineWorkspace

generateToyDataset: generateToyDataset.o libUtils.a
	$(LD) $(LDFLAGS) generateToyDataset.o $(LIBS) $(EXPLLINKLIBS) -o generateToyDataset

interpolatePdf: interpolatePdf.o libUtils.a
	$(LD) $(LDFLAGS) interpolatePdf.o $(LIBS) $(EXPLLINKLIBS) -o interpolatePdf

fritSignal: fritSignal.o libUtils.a
	$(LD) $(LDFLAGS) fritSignal.o $(LIBS) $(EXPLLINKLIBS) -o fritSignal

higgsSignalToPDF: higgsSignalToPDF.o libUtils.a
	$(LD) $(LDFLAGS) higgsSignalToPDF.o $(LIBS) $(EXPLLINKLIBS) -o higgsSignalToPDF

computeEff: computeEff.o libUtils.a
	$(LD) $(LDFLAGS) computeEff.o $(LIBS) $(EXPLLINKLIBS) -o computeEff

computeCrystalBallSyst: computeCrystalBallSyst.o  libUtils.a
	$(LD) $(LDFLAGS) computeCrystalBallSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeCrystalBallSyst

computeKeysPdfSyst: computeKeysPdfSyst.o  libUtils.a
	$(LD) $(LDFLAGS) computeKeysPdfSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeKeysPdfSyst

generateKeysPdfWorkspaces: generateKeysPdfWorkspaces.o  libUtils.a
	$(LD) $(LDFLAGS) generateKeysPdfWorkspaces.o $(LIBS) $(EXPLLINKLIBS) -o generateKeysPdfWorkspaces

computeJECSyst: computeJECSyst.o  libUtils.a
	$(LD) $(LDFLAGS) computeJECSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeJECSyst

computeJECSystWithToys: computeJECSystWithToys.o  libUtils.a
	$(LD) $(LDFLAGS) computeJECSystWithToys.o $(LIBS) $(EXPLLINKLIBS) -o computeJECSystWithToys

computeBkgSyst: computeBkgSyst.o  libUtils.a
	$(LD) $(LDFLAGS) computeBkgSyst.o $(LIBS) $(EXPLLINKLIBS) -o computeBkgSyst

computeSigmaRef: computeSigmaRef.o libUtils.a
	$(LD) $(LDFLAGS) computeSigmaRef.o $(LIBS) $(EXPLLINKLIBS) -o computeSigmaRef

runLikelihoodScan: runLikelihoodScan.o 
	$(LD) $(LDFLAGS) $^ $(DYNAMIC_LIBS) $(EXPLLINKLIBS) -o runLikelihoodScan

treatToyStuff: treatToyStuff.o  libUtils.a
	$(LD) $(LDFLAGS) treatToyStuff.o $(LIBS) $(EXPLLINKLIBS) -o treatToyStuff

extractLimits: extractLimits.o  libUtils.a
	$(LD) $(LDFLAGS) extractLimits.o $(LIBS) $(EXPLLINKLIBS) -o extractLimits

drawLimitCurve: drawLimitCurve.o  libUtils.a
	$(LD) $(LDFLAGS) drawLimitCurve.o $(LIBS) $(EXPLLINKLIBS) -o drawLimitCurve

libUtils.a: Utils.o
	@ar rcs $@ $^

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

# Clean everythingg
clean:
	@rm -f core* ${DEPS} ${OBJS} lib/*.o lib/*.d
	@rm -f *.a

# If we aren't doing a 'gmake clean' then make the dependencies
ifneq ($(MAKECMDGOALS), clean)
-include $(DEPS)
endif

