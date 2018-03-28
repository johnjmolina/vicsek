#ifndef COMMON_HXX
#define COMMON_HXX

#include "Array.h"
#include "Utils.hxx"

typedef unsigned int uint;
typedef Array::array1<int>::opt    vec_i;
typedef Array::array2<int>         nvec_i;

typedef Array::array1<double>::opt vec;
typedef Array::array2<double>      nvec;

typedef Array::array1<float>::opt  vec_f;
typedef Array::array2<float>       nvec_f;

// distance & PBC functions
using MDutils::distance2D;
using MDutils::pbc2D;
using Constants::PI;
using std::string;
using std::cerr;
using std::endl;
using std::exception;

// System flags
namespace NOISE{
  enum NOISE{
    ANGULAR,
    VECTORIAL,
    NONE,
    _SIZE_ 
  };
}
namespace OUTPUT{
  enum OUTPUT{
    SILO,
    GSD,
    HDF5,
    _SIZE_
  };
}
namespace LINK{
  enum LINK{
    OFF,
    ON,
    AUTO,
    MANUAL,
    SORTED,
    _SIZE_
  };
}
namespace INIT{
  enum INIT{
    RANDOM,
    RESTART,
    _SIZE_
  };
}

template<> struct enum_traits<NOISE::NOISE>: enum_traiter<NOISE::NOISE, NOISE::_SIZE_>{};
template<> struct enum_traits<LINK::LINK>  : enum_traiter<LINK::LINK,    LINK::_SIZE_> {};
template<> struct enum_traits<INIT::INIT>  : enum_traiter<INIT::INIT,    INIT::_SIZE_> {};
template<> struct enum_traits<OUTPUT::OUTPUT> : enum_traiter<OUTPUT::OUTPUT, OUTPUT::_SIZE_>{};

// System parameters
const int    DIM=2;            // system dimensions
const int    BUFFER_SIZE = 128;

// Simulation parameters
extern enum_set<NOISE::NOISE>   SW_NOISE;  // type of noise
extern enum_set<LINK::LINK>     SW_LINK;   // with link list?
extern enum_set<INIT::INIT>     SW_INIT;   // initialize
extern enum_set<OUTPUT::OUTPUT> SW_OUTPUT; // output format

extern int    NUM_SNAP;                    // number of snapshots
extern int    NUM_FRAMES_SNAP;             // number of frames / snapshot
extern int    MAX_NUM_STEPS;               // = num_dump * num_frames
extern int    ts;                          // time step


// Buckingham PI (adimensional) parameters
extern std::vector<double> PI_NOISE;
extern std::vector<double> PI_BOX;
extern double PI_VEL;
extern double PI_RHO;

extern double DT;                     // unit of time
extern double R0;                     // unit of length

extern double ETA_ANG;                // temperature angular noise
extern double ETA_VEC;                // temperature vectorial noise
extern double RHO;                    // density
extern double U;                      // particle speed
extern double UDT;                    // particle speed * time step

extern int    NUMP;                   // number of particles
extern double BOXL[DIM];              // box dimension

extern double RISQ;                   // alignment radius

// particle properties
extern nvec   pos;                    // particle position
extern nvec   vel;                    // particle orientation (unit velocity vector)

extern vec    theta;                  // particle orientation
extern vec    xi;                     // particle angle fluctuations
extern nvec   vxi;                    // particle orientation fluctuations
extern nvec   vlop;                   // vector local order parameter
extern vec    lop;                    // local order parameter
extern vec    sum_sin;                // average sin theta over neighbors
extern vec    sum_cos;                // average cos theta over neighbors
extern vec_i  sum_cnt;                // number of neighbors
#endif
