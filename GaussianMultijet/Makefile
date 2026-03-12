RK_LIBDIR = /usr/local/lib

OFILES = stringUtils.o genClusterAssignments.o permutation.o mathUtils.o

PROGRAMS = testClusteringSequence.cc testMultiParticleCollection.cc \
           testEventProcessingFramework.cc gaussianMultijetStudy.cc \
	   gaussianMultijetClustComp.cc pythiaJetStudy.cc

OPTIMIZE = -O3
# OPTIMIZE = -O0 -ggdb
INCLUDES = -I. -I$(RK_LIBDIR)/../include
# ROOTLIBS=`root-config --glibs`
# ROOTLIBDIR = `root-config --libdir`
CPPFLAGS = $(OPTIMIZE) $(INCLUDES) -Wall -W -Werror

LIBS = -L$(RK_LIBDIR) -lsjet -lcnpy -lrk -lpythia8 -lfastjet $(ROOTLIBS) -ldl -lm -lz

%.o : %.cc
	g++ -c $(CPPFLAGS) -fPIC -MD $< -o $@
	@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

BINARIES = $(PROGRAMS:.cc=)

all: $(BINARIES)

$(BINARIES): % : %.o $(OFILES); g++ $(OPTIMIZE) -fPIC -o $@ $^ $(LIBS) --std=c++11

clean:
	rm -f $(BINARIES) core core.* *.o *.d *~

-include $(OFILES:.o=.d)
-include $(PROGRAMS:.cc=.d)
