#include "Common.hxx"

// Simulation parameters
enum_set<NOISE::NOISE>   SW_NOISE;  // type of noise
enum_set<LINK::LINK>     SW_LINK;   // with link list?
enum_set<INIT::INIT>     SW_INIT;   // initialize
enum_set<OUTPUT::OUTPUT> SW_OUTPUT; // output format

int    NUM_SNAP;                    // number of snapshots
int    NUM_FRAMES_SNAP;             // number of frames / snapshot
int    MAX_NUM_STEPS;               // = num_dump * num_frames
int    ts;                          // time step


// Buckingham PI (adimensional) parameters
std::vector<double> PI_NOISE;
std::vector<double> PI_BOX;
double              PI_VEL;
double              PI_RHO;

double DT;                     // unit of time
double R0;                     // unit of length

double ETA_ANG;                // temperature angular noise
double ETA_VEC;                // temperature vectorial noise
double RHO;                    // density
double U;                      // particle speed
double UDT;                    // particle speed * time step

int    NUMP;                   // number of particles
double BOXL[DIM];              // box dimensions

double RISQ;                   // alignment radius

// particle properties
nvec   pos;                    // particle position
nvec   vel;                    // particle orientation (unit velocity vector)

vec    theta;                  // particle orientation
vec    xi;                     // particle angle fluctuations
nvec   vxi;                    // particle orientation fluctuations
nvec   vlop;                   // vector local order parameter
vec    lop;                    // 
vec    sum_sin;                // average sin theta over neighbors
vec    sum_cos;                // average cos theta over neighbors
vec_i  sum_cnt;                // number of neighbors
