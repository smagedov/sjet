# Change this variable to where you installed pythia, sjet, rk
LOCALLIBS = /usr/local/lib

FASTJETCONFIG   := fastjet-config

FASTJETPREFIX   := $(shell $(FASTJETCONFIG) --prefix)
FASTJETLIBS     := $(shell $(FASTJETCONFIG) --libs)
FASTJETCPPFLAGS := $(shell $(FASTJETCONFIG) --cxxflags)

OPTIMIZE = -O3
CPPFLAGS = $(OPTIMIZE) -I. $(FASTJETCPPFLAGS) -I$(LOCALLIBS)/../include
LIBS = -L$(LOCALLIBS) -lpythia8 $(FASTJETLIBS) -lsjet -lrk -ldl -lm

PROGRAMS = slow_example.cc slow_pythia.cc fastjet_example.cc fastjet_pythia.cc

%.o : %.cc
	g++ -c $(CPPFLAGS) -fPIC -MD $< -o $@
	@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

BINARIES = $(PROGRAMS:.cc=)

all: $(BINARIES)

$(BINARIES): % : %.o $(OFILES); g++ $(OPTIMIZE) -fPIC -o $@ $^ $(LIBS)

clean:
	rm -f *.so *.o *.d *~ $(BINARIES) core.*

-include $(OFILES:.o=.d)
-include $(PROGRAMS:.cc=.d)
