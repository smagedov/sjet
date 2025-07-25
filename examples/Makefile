FASTJETCONFIG   := fastjet-config

FASTJETPREFIX   := $(shell $(FASTJETCONFIG) --prefix)
FASTJETLIBS     := $(shell $(FASTJETCONFIG) --libs)
FASTJETLIBDIR   := $(FASTJETPREFIX)/lib
FASTJETCPPFLAGS := $(shell $(FASTJETCONFIG) --cxxflags)

LOCALLIBS = /home/smagedov/local/lib
OPTIMIZE = -O2
CPPFLAGS = $(OPTIMIZE) -I. $(FASTJETCPPFLAGS) -I../include -I$(LOCALLIBS)/../include -I/usr/local/root/include
LIBS = -L$(LOCALLIBS) -lrk -lm -L../lib -lpythia8 -ldl -L/home/smagedov/local/lib -lfastjet -L/usr/local/root/lib -lCore -L/usr/local/root/lib -lGui -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lROOTNTuple -lMultiProc -lROOTDataFrame -lROOTNTupleUtil -ldl

OFILES = DiffusionDistFastJet.o diffusionFactor2.o

PROGRAMS = slow_example.cc slow_pythia.cc fastjet_example.cc fastjet_pythia.cc

%.o : %.cc
	g++ -c $(CPPFLAGS) -fPIC -MD $< -o $@
	@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

BINARIES = $(PROGRAMS:.cc=)

all: libfastjetdiffusion.so $(BINARIES)

$(BINARIES): % : %.o $(OFILES); g++ $(OPTIMIZE) -fPIC -o $@ $^ -w -L. -lfastjetdiffusion $(LIBS) -std=c++11 -pedantic -W -Wall -Wshadow -pthread -Wl,-rpath,../li -Wl,-rpath,/home/smagedov/local/lib -Wl,-rpath,/usr/local/root/lib -pthread -std=c++17 -m64 -Wl,-rpath,/usr/local/root/lib -pthread -rdynamic -DPY8ROOT

libfastjetdiffusion.so: $(OFILES)
	g++ -shared -o $@ $^ $(LIBS)

clean:
	rm -f *.so *.o *.d *~ $(BINARIES) core.*

-include $(OFILES:.o=.d)
-include $(PROGRAMS:.cc=.d)
