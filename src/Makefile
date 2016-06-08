# A Generic Makefile for compiling ROOT programs
# R. Michaels, rom@jlab.org, Aug 2001  See also README !!
# Version of this release
#
bin	      = /w/hallb-scifs1a/primex/yangz/work/jobs/bin/
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
CXXFLAGS      = -O -m64 -std=c++0x -Wall -frtti -fPIC \
		-lMinuit -DLINUXVERS -I$(ROOTSYS)/include/ -I../inc



# Linux with g++
INCLUDES      = -I$(ROOTSYS)/include
CXX           = g++
LD            = g++
LDFLAGS       =

LIBS          = $(ROOTLIBS) -lTreePlayer
GLIBS         = $(ROOTGLIBS) -lTreePlayer -L/opt/X11/lib -lXpm -lX11

ALL_LIBS =  $(GLIBS) $(LIBS)
#ALL_LIBS =  $(LIBS)

$(O): $(O).o primex.o
	$(CXX) $(CXXFLAGS) -o $(bin)/$@ $(O).o primex.o $(ALL_LIBS)

all: rotatedmass rotatedmass_mc rotatedmass_omg fast_fitmodel_mc fast_fitmodel_mc_correct altfit alt_model_correct plotalt dndt_fit btcorr

combine: combine.o
	$(CXX) $(CXXFLAGS) combine.o $(ALL_LIBS) -o $(bin)/$@

change_name: change_name.o primex.o
	$(CXX) $(CXXFLAGS) change_name.o $(ALL_LIBS) -o $(bin)/$@

rotatedmass: rotatedmass.o primex.o
	$(CXX) $(CXXFLAGS) rotatedmass.o primex.o $(ALL_LIBS) -o $(bin)/$@  

rotatedmass_mc: rotatedmass_mc.o primex.o
	$(CXX) $(CXXFLAGS) rotatedmass_mc.o primex.o $(ALL_LIBS) -o $(bin)/$@  

rootperf: rootperf.o primex.o
	$(CXX) $(CXXFLAGS) rootperf.o primex.o $(ALL_LIBS) -o $(bin)/$@  

rotatedmass_omg: rotatedmass_omg.o primex.o
	$(CXX) $(CXXFLAGS) rotatedmass_omg.o primex.o $(ALL_LIBS) -o $(bin)/$@  

fast_fitmodel_mc: fast_fitmodel_mc.o primex.o
	$(CXX) $(CXXFLAGS) fast_fitmodel_mc.o primex.o $(ALL_LIBS) -o $(bin)/$@  

fast_fitmodel_mc_correct: fast_fitmodel_mc_correct.o primex.o alt_model_fit.o
	$(CXX) $(CXXFLAGS) fast_fitmodel_mc_correct.o primex.o alt_model_fit.o $(ALL_LIBS) -o $(bin)/$@  

altfit: altfit.o primex.o alt_model_fit.o
	$(CXX) $(CXXFLAGS) altfit.o primex.o alt_model_fit.o $(ALL_LIBS) -o $(bin)/$@  

altfit.0.8: altfit.o primex.o alt_model_fit_0.8.o
	$(CXX) $(CXXFLAGS) altfit.o primex.o alt_model_fit_0.8.o $(ALL_LIBS) -o $(bin)/$@  

altfit.1.2: altfit.o primex.o alt_model_fit_1.2.o
	$(CXX) $(CXXFLAGS) altfit.o primex.o alt_model_fit_1.2.o $(ALL_LIBS) -o $(bin)/$@  

alt_model_correct: alt_model_correct.o primex.o alt_model_fit.o
	$(CXX) $(CXXFLAGS) alt_model_correct.o primex.o alt_model_fit.o $(ALL_LIBS) -o $(bin)/$@  

alt_model_correct_0.8: alt_model_correct.o primex.o alt_model_fit_0.8.o
	$(CXX) $(CXXFLAGS) alt_model_correct.o primex.o alt_model_fit_0.8.o $(ALL_LIBS) -o $(bin)/$@  

alt_model_correct_1.2: alt_model_correct.o primex.o alt_model_fit_1.2.o
	$(CXX) $(CXXFLAGS) alt_model_correct.o primex.o alt_model_fit_1.2.o $(ALL_LIBS) -o $(bin)/$@  

plotalt: plotalt.o primex.o
	$(CXX) $(CXXFLAGS) plotalt.o primex.o $(ALL_LIBS) -o $(bin)/$@  

dndt_fit: dndt_fit.o primex.o
	$(CXX) $(CXXFLAGS) dndt_fit.o primex.o $(ALL_LIBS) -o $(bin)/$@  

btcorr: btcorr.o primex.o alt_model_fit.o
	$(CXX) $(CXXFLAGS) -o $(bin)/$@ btcorr.o primex.o alt_model_fit.o $(ALL_LIBS)

test: test.o primex.o
	$(CXX) $(CXXFLAGS) -o $(bin)/$@ test.o primex.o $(ALL_LIBS)

rotatedmass.o: rotatedmass.cpp
	$(CXX) $(CXXFLAGS) -c rotatedmass.cpp 

rotatedmass_mc.o: rotatedmass_mc.cpp
	$(CXX) $(CXXFLAGS) -c rotatedmass_mc.cpp 

rootperf.o: rootperf.cpp
	$(CXX) $(CXXFLAGS) -c rootperf.cpp 

rotatedmass_omg.o: rotatedmass_omg.cpp
	$(CXX) $(CXXFLAGS) -c rotatedmass_omg.cpp 

fast_fitmodel_mc.o: fast_fitmodel_mc.cpp
	$(CXX) $(CXXFLAGS) -c fast_fitmodel_mc.cpp 

fast_fitmodel_mc_correct.o: fast_fitmodel_mc_correct.cpp
	$(CXX) $(CXXFLAGS) -c fast_fitmodel_mc_correct.cpp 

altfit.o: altfit.cpp
	$(CXX) $(CXXFLAGS) -c altfit.cpp 

alt_model_correct.o: alt_model_correct.cpp
	$(CXX) $(CXXFLAGS) -c alt_model_correct.cpp 

primex.o: primex.cpp
	$(CXX) $(CXXFLAGS) -c primex.cpp 

dndt_fit.o: dndt_fit.cpp
	$(CXX) $(CXXFLAGS) -c dndt_fit.cpp 

alt_model_fit.o: alt_model_fit.cpp
	$(CXX) $(CXXFLAGS) -c alt_model_fit.cpp 

alt_model_fit_0.8.o: alt_model_fit_0.8.cpp
	$(CXX) $(CXXFLAGS) -c alt_model_fit_0.8.cpp 

alt_model_fit_1.2.o: alt_model_fit_1.2.cpp
	$(CXX) $(CXXFLAGS) -c alt_model_fit_1.2.cpp 

plotalt.o: plotalt.cpp
	$(CXX) $(CXXFLAGS) -c plotalt.cpp 

btcorr.o: btcorr.cpp
	$(CXX) $(CXXFLAGS) -c btcorr.cpp

clean:
	rm -f *.o *.d $(bin)/*

###
.SUFFIXES:
.SUFFIXES: .c .cc .cpp .C .o .d

%.o:    %.cpp
	$(CXX) $(CXXFLAGS) -c $<

%.d:    %.cpp
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(CXX) -MM $(CXXFLAGS) -c $< \
                | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@; \
                [ -s $@ ] || rm -f $@'

-include $(DEPS)
