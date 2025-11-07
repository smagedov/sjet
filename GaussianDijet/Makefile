#RK_LIBDIR = /usr/local/lib
LOCALLIBS = /uscms/home/wejiejin/work/rk_install/lib
PYTHIALIB = /uscms/home/wejiejin/work/pythia8315/lib
PYTHIAINC = /uscms/home/wejiejin/work/pythia8315/include
FASTJETLIB = /uscms/home/wejiejin/work/pythia8315/lib
OFILES = stringUtils.o genClusterAssignments.o permutation.o mathUtils.o

PROGRAMS = testClusteringSequence.cc testMultiParticleCollection.cc \
           testEventProcessingFramework.cc gaussianDijetStudy.cc \
           gaussianJetStudy.cc pythiaTTbarStudy.cc

OPTIMIZE = -O2
# OPTIMIZE = -O0 -ggdb


CPPFLAGS = $(OPTIMIZE) -I. -I../include -I$(LOCALLIBS)/../include -I$(PYTHIAINC) -I`root-config --incdir`
ROOTLIBS=`root-config --glibs`
ROOTLIBDIR = `root-config --libdir`
LIBS = -L$(LOCALLIBS) -lrk -lm -L$(PYTHIALIB) -lpythia8 -ldl -L$(FASTJETLIB) -lfastjet $(ROOTLIBS) -lsjet




#INCLUDES = -I. -I$(RK_LIBDIR)/../include
#CPPFLAGS = $(OPTIMIZE) $(INCLUDES) -Wall -W -Werror

#LIBS = -L$(RK_LIBDIR) -lsjet -lrk -lm

%.o : %.cc
	g++ -c $(CPPFLAGS) -fPIC -MD $< -o $@
	@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

BINARIES = $(PROGRAMS:.cc=)

all: $(BINARIES)

#$(BINARIES): % : %.o $(OFILES); g++ $(OPTIMIZE) -fPIC -o $@ $^ $(LIBS)


$(BINARIES): % : %.o $(OFILES); g++ $(OPTIMIZE) -fPIC -o $@ $^ -w $(LIBS) -std=c++11 -pedantic -W -Wall -Wshadow -pthread -Wl,-rpath,../li -Wl,-rpath,$(LOCALLIBS) -Wl,-rpath,$(ROOTLIBDIR) -Wl,-rpath,$(PYTHIALIB) -Wl,-rpath,$(FASTJETLIB) -pthread -std=c++17 -m64 -Wl,-rpath,$(ROOTLIBDIR) -pthread -rdynamic -DPY8ROOT

clean:
	rm -f $(BINARIES) core core.* *.o *.d *~

-include $(OFILES:.o=.d)
-include $(PROGRAMS:.cc=.d)
