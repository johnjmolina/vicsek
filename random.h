#ifndef MACRO_HXX
#define MACRO_HXX

#include <cstdlib>
#include <iostream>
#include "Utils.hxx"

// random number generator
#include <ctime>
#include "mt19937ar.h"

// random number seeds
namespace  RND_SEED{
  enum RND_SEED{GIVEN, RANDOM};
}

// random number help functions
// initialize random number seed (if necessary)
inline void RA_SEED(const RND_SEED::RND_SEED &SW_SEED, const unsigned long &seed){
  using std::cerr;
  using std::endl;
  cerr << "# Using Marseene Twister 19937 generator\n";
  if(SW_SEED == RND_SEED::RANDOM){
    init_genrand(time(NULL));
    cerr << "# Random seed = time(NULL)\n";
  }else if(SW_SEED == RND_SEED::GIVEN){
    init_genrand(seed);
    cerr << "# Random seed = " << seed << endl;
  }else{
    cerr << "# Invalid seed\n";
    exit(1);
  }
}


//uniform random [0, x]
inline double RA_CC(const double &x){
  return (genrand_real1()*x);
}
//uniform random [a, b]
inline double RA_CC(const double &a, const double &b){
  return (ABS(b-a)*(genrand_real1()) + MIN(a, b));
}

//uniform random [0, x)
inline double RA_CO(const double &x){
  return (genrand_real2()*x);
}
//uniform random [a, b)
inline double RA_CO(const double &a, const double &b){
  return (ABS(b-a)*(genrand_real2()) + MIN(a, b));
}

//uniform random (0, x)
inline double RA_OO(const double &x){
  return (genrand_real3()*x);
}
//uniform random (a, b)
inline double RA_OO(const double &a, const double &b){
  return (ABS(b-a)*(genrand_real3()) + MIN(a, b));
}


#endif
