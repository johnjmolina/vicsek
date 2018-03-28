#ifndef OUTPUT_HXX
#define OUTPUT_HXX

#include <silo.h>
#include "SimIO.hxx"
#include "Common.hxx"
#include "gsd.h"



// output parameters
extern string dir_name;
extern string prj_name;

extern SimIO::SimIO             h5out;   // hdf5 output writer
extern SimIO::io_space h5scalar_space;  // dataspace for 1d data
extern SimIO::io_space h5vector_space;  // dataspace for 2d data

extern Array::array1<float>         float_buffer;
extern Array::array1<unsigned int>  uint_buffer;


void hdf5_init();
void hdf5_output(int frame_id);
void hdf5_end();

void silo_init();
void silo_output(int frame_id);
void silo_end();

void gsd_init();
void gsd_output(int frame_id);
void gsd_end();

#endif
