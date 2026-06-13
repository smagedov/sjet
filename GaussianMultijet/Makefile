RK_LIBDIR = /usr/local/lib

TORCH_DIR = /home/smagedov/libtorch

OFILES = stringUtils.o genClusterAssignments.o permutation.o mathUtils.o

PROGRAMS = testClusteringSequence.cc testMultiParticleCollection.cc \
           testEventProcessingFramework.cc gaussianMultijetStudy.cc \
	   gaussianMultijetClustComp.cc pythiaJetStudy.cc

OPTIMIZE = -O1 -g
# OPTIMIZE = -O0 -ggdb
INCLUDES = -I. -I$(RK_LIBDIR)/../include -I$(TORCH_DIR)/include -I$(TORCH_DIR)/include/torch/csrc/api/include
# ROOTLIBS=`root-config --glibs`
# ROOTLIBDIR = `root-config --libdir`
CPPFLAGS = $(OPTIMIZE) $(INCLUDES) #-Wall -W -Werror

LIBS = -L$(RK_LIBDIR) -lsjet -lcnpy -lrk -lpythia8 \
       -lfastjet $(ROOTLIBS) -ldl -lm -lz \
       -L$(TORCH_DIR)/lib -ltorch -ltorch_cpu -lc10

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
