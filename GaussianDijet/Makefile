RK_LIBDIR = /usr/local/lib

OFILES = stringUtils.o genClusterAssignments.o permutation.o mathUtils.o

PROGRAMS = testClusteringSequence.cc testMultiParticleCollection.cc \
           testEventProcessingFramework.cc gaussianDijetStudy.cc \
           gaussianJetStudy.cc

OPTIMIZE = -O3
# OPTIMIZE = -O0 -ggdb
INCLUDES = -I. -I$(RK_LIBDIR)/../include
CPPFLAGS = $(OPTIMIZE) $(INCLUDES) -Wall -W -Werror

LIBS = -L$(RK_LIBDIR) -lsjet -lrk -lm

%.o : %.cc
	g++ -c $(CPPFLAGS) -fPIC -MD $< -o $@
	@sed -i 's,\($*\.o\)[:]*\(.*\),$@: $$\(wildcard\2\)\n\1:\2,g' $*.d

BINARIES = $(PROGRAMS:.cc=)

all: $(BINARIES)

$(BINARIES): % : %.o $(OFILES); g++ $(OPTIMIZE) -fPIC -o $@ $^ $(LIBS)

clean:
	rm -f $(BINARIES) core core.* *.o *.d *~

-include $(OFILES:.o=.d)
-include $(PROGRAMS:.cc=.d)
