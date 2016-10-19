INC_ADETL= -I/opt/FuRSST_SOFTWARE/ADETL
INC_MESH= -I/opt/FuRSST_SOFTWARE/DiscreteDomain
INC_TAB=-I/opt/FuRSST_SOFTWARE/LookupTables
INC_SOL=-I/opt/FuRSST_SOFTWARE/Solvers
INC_MKL=-I/opt/intel/mkl/include
INC= -I./src $(INC_ADETL) $(INC_MESH) $(INC_TAB) $(INC_SOL) $(INC_MKL) $(INC_SAMG)
SAMGLIBDIR=/opt/samg27a2_omp-linuxX64-lstdc++6-ifort12.1.3.293
IFORT_RTS-LIBS=/opt/intel/lib/intel64/libiomp5.so
LIBAMG=$(SAMGLIBDIR)/libamg.so $(IFORT_RTS-LIBS)
AMG_DEFS=-DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE
INC_SAMG=-I$(SAMGLIBDIR)

MKLRTDIR = /opt/intel/mkl/lib/intel64
LIBS= $(LIBAMG) -L$(MKLRTDIR) -Wl,--start-group $(MKLRTDIR)/libmkl_intel_lp64.a $(MKLRTDIR)/libmkl_core.a $(MKLRTDIR)/libmkl_intel_thread.a -Wl,--end-group -lpthread -lm

CFLAGS= -g -m64 -Ofast -std=c++11 -openmp $(AMG_DEFS)
SRC= ./src/*.cpp 

all:
	icpc $(CFLAGS) $(SRC) $(INC) -o ./run $(LIBS)

clean:
	rm -rf ./run ./src/*~ ./*~ ./output/*~ ./input/*~

