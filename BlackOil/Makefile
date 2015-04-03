INC_ADETL= -I../../../Library/ADETL
#INC_SUPERLU= -I/opt/SuperLU_SERIAL/include/
#INC= -I./src $(INC_ADETL) $(INC_SUPERLU)
MKLROOT = /opt/intel/mkl/
#PARDISO = -I/opt/intel/mkl/include/
STATIC= -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm

DYNAMIC= -Wl,--no-as-needed -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm

OPTION= -m64 -I$(MKLROOT)/include
#LIBDIR= -L/opt/SuperLU_SERIAL/lib



CFLAGS= -Wall -g -O3 -std=c++11
SRC= *.cpp 

all:
	g++ $(CFLAGS) $(INC_ADETL) $(OPTION) $(SRC) -o ./run $(DYNAMIC)

clean:
	rm -rf ./run ./*~
