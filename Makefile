ifeq ($(ENV), CLANG)
     CC       = clang
     CXX      = clang++
     COPT     = -O3 -I/usr/local/include
     CCOPT    = $(COPT) -std=c++11 
     LINKS    = -lm -L/usr/local/lib -lsiloh5 -lhdf5_hl -lhdf5
     EXE      = sopt
endif
ifeq ($(ENV), GCC)
     SILO_DIR = /opt/silo
     HDF5_DIR = /opt/hdf5
     CC       = gcc-7
     CXX      = g++-7
     COPT     = -DNDEBUG -O3 -I$(HDF5_DIR)/include -I$(SILO_DIR)/include -I/opt/local/include -I/usr/local/include
#     COPT     = -DNDEBUG -O3 -I/opt/hdf5/include
     CCOPT    = $(COPT) -std=c++11 
     LINKS    = -lm -L$(SILO_DIR)/lib -L$(HDF5_DIR)/lib -lsiloh5 -lhdf5_hl -lhdf5
     EXE      = sopt
endif
ifeq ($(ENV), GCC_OMP)
     SILO_DIR = /opt/silo
     HDF5_DIR = /opt/hdf5
     CC       = gcc-6
     CXX      = g++-6
     COPT     = -DNDEBUG -O3 -I$(HDF5_DIR)/include -I$(SILO_DIR)/include -fopenmp -I/opt/local/include -I/usr/local/include
     CCOPT    = $(COPT) -std=c++11 
     LINKS    = -lm -L$(SILO_DIR)/lib -L$(HDF5_DIR)/lib -lsiloh5 -lhdf5_hl -lhdf5
     EXE      = popt
endif
ifeq ($(ENV), ICC)
     CC       = icc
     CXX      = icpc

     HDF5DIR  = /opt/hdf5.1.8
     BOOSTDIR = /usr/local/boost_1_53_0
     COPT     = -DNDEBUG -O3 -I/usr/local/include -I$(HDF5DIR)/include -I$(BOOSTDIR) \
	-w0 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2
     CCOPT    = $(COPT) -std=c++11
     LINKS    = -L/usr/local/lib -L$(HDF5DIR)/lib -lhdf5_hl -lhdf5 -lm -lstdc++
     EXE      = sopt
endif
ifeq ($(ENV), ICC_OMP)
     MKL_DIR  = /home/opt/intel/compilers_and_libraries/linux
     MKL_LIB_PATH = $(MKL_DIR)/lib/intel64
     MKL_INCLUDE_PATH = $(MKL_DIR)/include

     HDF5DIR  = /opt/hdf5.1.8
     BOOSTDIR = /usr/local/boost_1_53_0

     CC       = icc
     CXX      = icpc

     COPT     = -DNDEBUG -O3 -ip -qopenmp -parallel -w0 -xSSSE3 -axAVX,SSE4.2,SSE4.1,SSSE3,SSE3,SSE2 \
		-I/usr/local/include -I$(HDF5DIR)/include -I$(BOOSTDIR)
     CCOPT    = $(COPT)  -std=c++11 -I$(MKL_INCLUDE_PATH)
     LINKS    = -lstdc++ -L/usr/local -L$(HDF5DIR)/lib -lhdf5_hl -lhdf5 \
		-L$(MKL_LIB_PATH) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm
     EXE      = popt
endif
CFLAGS   = $(CCOPT)

OBJS =  Utils.o \
	Common.o \
	Output.o \
	Vicsek.o \
	gsd.o \
	mt19937ar.o

## Implicit rules

.SUFFIXES: .c .cxx .o .out

VICSEK = vicsek

## Build rules
all: $(VICSEK) 

$(VICSEK): $(OBJS)
	$(CXX) $(OBJS) -o $(VICSEK).$(EXE) $(CFLAGS) $(LINKS)

## Compile

.cxx.o: 
	$(CXX) -c $< $(CCOPT) -o $@

.c.o: 
	$(CC) -c $< $(COPT) -o $@

## Clean
clean:
	rm -f *.o *.dat *~
cleanall:
	rm -f *.o *.dat *.h5 *~ *.out *.err *.x *.sopt *.popt *.xdmf

