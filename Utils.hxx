#ifndef UTILS_HPP
#define UTILS_HPP

#include <cmath>
#include <limits>
#include <bitset>
#include <stdexcept>
#include <cassert>
#include <vector>
#include <cctype>
#include <string>
#include <iostream>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>


 
// openmp support
#ifdef _OPENMP
#include <omp.h>
#else
#ifndef __clang__
#include <chrono>
#else
#include <ctime>
#endif
#endif

namespace Constants{
  static const double PI    = M_PI;
  static const double TWOPI = 2.0*M_PI;
  static const double PI_HALF = PI / 2.0;
  static const double SIXTH_ROOT_TWO = pow(2.0, 1.0/6.0);
  static const double THIRD_ROOT_TWO = SIXTH_ROOT_TWO*SIXTH_ROOT_TWO;
  static const double EPSILON_MP= std::numeric_limits<double>::epsilon(); //1.0e-15
  static const double MAX_MP    = std::numeric_limits<double>::max();
  static const double TOL_MP      =10.0*EPSILON_MP;   //1.0e-14
  static const double BIG_TOL_MP  =1.0e3*TOL_MP;      //1.0e-11
  static const double BIGG_TOL_MP =1.0e3*BIG_TOL_MP;  //1.0e-8
  static const double HUGE_TOL_MP =1.0e3*BIGG_TOL_MP; //1.0e-5
}
inline void init_threads(){
  using std::cerr;
  using std::endl;
#ifdef _OPENMP
  int nthreads, tid, procs, maxt, inpar, dynamic, nested;
#pragma omp parallel private(nthreads, tid)
  {
    tid = omp_get_thread_num();
    if(tid == 0){
      procs   = omp_get_thread_num();
      nthreads= omp_get_num_threads();
      maxt    = omp_get_max_threads();      
      inpar   = omp_in_parallel();      
      dynamic = omp_get_dynamic();
      nested  = omp_get_nested();

      cerr << "# " << endl;
      cerr << "# OMP RUNTIME :" << endl;
      cerr << "# Number of processors        = " << procs    << endl;
      cerr << "# Number of threads           = " << nthreads << endl;
      cerr << "# Max threads                 = " << maxt     << endl;
      cerr << "# In parallel ?               = " << inpar    << endl;
      cerr << "# Dynamic thread enabled ?    = " << dynamic  << endl;
      cerr << "# Nested parallelism enabled? = " << nested   << endl;
      cerr << "# " << endl;
    }
  }
#endif
}
inline std::string tolower(const std::string &dmy){
  std::string low(dmy);
  for(unsigned int i = 0; i < low.size(); i++) low[i] = tolower(low[i]);
  return low;
}
inline std::string toupper(const std::string &dmy){
  std::string upp(dmy);
  for(unsigned int i = 0; i < upp.size(); i++) upp[i] = toupper(upp[i]);
  return upp;
}
inline void run_assert(const bool& p){
  if(!p) throw std::runtime_error("runtime assert failed!");
}
inline void run_assert(const bool& p, const std::string& arg){
  if(!p) throw std::runtime_error(arg);
}
inline void run_assert(const bool& p, const char* arg){
  if(!p) throw std::runtime_error(arg);
}

template<typename T>
inline bool non_zero_mp(const T &a){
  return (a > Constants::TOL_MP) || (-Constants::TOL_MP > a);
}
template <typename T>
inline bool zero_mp(const T &a){
  return (a <= Constants::TOL_MP && (-Constants::TOL_MP <= a));
}
template <typename T>
inline bool positive_mp(const T &a){
  return a > Constants::TOL_MP;
}
template <typename T>
inline bool negative_mp(const T &a){
  return -Constants::TOL_MP > a;
}

// nearest integer
template<typename T>
inline int Nint(const T &a){
  return (a >= 0.0 ? static_cast<int>(a + 0.5) : static_cast<int>(a - 0.5));
}

// square
template<typename T>
inline T SQ(const T &a) {
  return a*a;
}

template<typename T>
inline T POW3(const T &a){
  return a*a*a;
}

// minimum
template<typename T>
inline T MIN(const T &a, const T &b){
  return ( a <= b ? a : b);
}

// maximum
template<typename T>
inline T MAX(const T &a, const T &b){
  return ( a >= b ? a : b);
}

// absolute value
template<typename T>
inline T ABS(const T &a){
  return ( a >= 0 ? a : -a);
}

template<typename T>
inline void SWAP(T &a, T &b){
  T a_cp = a;
  a = b;
  b = a;
}

inline double equal_tol(const double &a, const double &b, const double &rtol, const double &atol = Constants::TOL_MP){
  if(a == b) return true;

  double diff = ABS(a - b);
  if(diff <= atol) return true;

  double eps = MAX(ABS(a), ABS(b))*rtol;
  return (diff <= eps ? true : false);
}
inline bool equal_mp(const double &a, const double &b){
  if(a == b) return true;
  double eps = (MAX(ABS(a), ABS(b)) + 10.0)*Constants::EPSILON_MP;  
  return (ABS(a - b) <= eps ? true : false);
}

inline bool dirCheckMake(const char* dname){
  DIR* dtest;
  if((dtest=opendir(dname))==NULL){
    char dmy_cmd[256];
    snprintf(dmy_cmd, 256, "mkdir %s", dname);
    system(dmy_cmd);
    DIR* dnew;
    if((dnew=opendir(dname)) == NULL){
      fprintf(stderr, "Error: failed to create directory\n"); 
      return false;
    }
    closedir(dnew);
  }else{
    closedir(dtest);
  }
  return true;
}

struct WallTimer{
#ifndef _OPENMP
#ifndef __clang__
  std::chrono::time_point<std::chrono::system_clock> t_start;
  inline void start(){t_start = std::chrono::system_clock::now();}
  inline double stop(){
    std::chrono::time_point<std::chrono::system_clock> t_end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = t_end - t_start;
    return elapsed.count();
  }
#else
  time_t t_start;  
  inline void start(){time(&t_start);}
  inline double stop(){
    time_t t_end;
    time(&t_end);
    return static_cast<double>(difftime(t_end, t_start));
  }
#endif
#else
  double t_start;
  inline void start(){t_start = omp_get_wtime();}
  inline double stop(){return (omp_get_wtime() - t_start);}
#endif
};

// MD utilities for orthorombic unit cell with origin at 0
namespace MDutils{
  // compute pbc distance
  inline double distance(double *r12,
			 const double *r1, const double *r2,
			 const double *lbox, const int &dim){
    double dmy_dist = 0.0;
    for(int d = 0; d < dim; d++){
      r12[d] = r2[d] - r1[d];
      r12[d] -= static_cast<double>(Nint(r12[d]/lbox[d])) * lbox[d];
      dmy_dist += SQ(r12[d]);
    }
    return dmy_dist;
  }
  inline double distance(double *r12,
			 const double *r1, const double *r2,
			 const double *lbox, const int *pbc_boundary, const int &dim){
    double dmy_dist = 0.0;
    for(int d = 0; d < dim; d++){
      r12[d] = r2[d] - r1[d];
      if(pbc_boundary[d]) r12[d] -= static_cast<double>(Nint(r12[d]/lbox[d])) * lbox[d];
      dmy_dist += SQ(r12[d]);
    }
    return dmy_dist;
  }
  inline double distance2D(double *r12, const double *lbox){
    r12[0] -= static_cast<double>(Nint(r12[0]/lbox[0])) * lbox[0];
    r12[1] -= static_cast<double>(Nint(r12[1]/lbox[1])) * lbox[1];
    return (SQ(r12[0]) + SQ(r12[1]));
  }
  inline double distance2D(double *r12, const double *r1, const double *r2, const double *lbox){
    r12[0] = r2[0] - r1[0];
    r12[0] -= static_cast<double>(Nint(r12[0]/lbox[0])) * lbox[0];    

    r12[1] = r2[1] - r1[1];
    r12[1] -= static_cast<double>(Nint(r12[1]/lbox[1])) * lbox[1];
    
    return (SQ(r12[0]) + SQ(r12[1]));
  }
  inline double distance2D(double *r12, const double *r1, const double *r2, const double *lbox,
			   const int *pbc_boundary){
    r12[0] = r2[0] - r1[0];
    r12[1] = r2[1] - r1[1];
    
    if(pbc_boundary[0]) r12[0] -= static_cast<double>(Nint(r12[0]/lbox[0])) * lbox[0];
    if(pbc_boundary[1]) r12[1] -= static_cast<double>(Nint(r12[1]/lbox[1])) * lbox[1];
    
    return (SQ(r12[0]) + SQ(r12[1]));
  }
  inline double distance3D(double *r12, const double *lbox){
    r12[0] -= static_cast<double>(Nint(r12[0]/lbox[0])) * lbox[0];
    r12[1] -= static_cast<double>(Nint(r12[1]/lbox[1])) * lbox[1];
    r12[2] -= static_cast<double>(Nint(r12[2]/lbox[2])) * lbox[2];
    return (SQ(r12[0]) + SQ(r12[1]) + SQ(r12[2]));
  }
  inline double distance3D(double *r12, const double *r1, const double *r2, const double *lbox){
    r12[0] = r2[0] - r1[0];
    r12[0] -= static_cast<double>(Nint(r12[0]/lbox[0])) * lbox[0];
    
    r12[1] = r2[1] - r1[1];
    r12[1] -= static_cast<double>(Nint(r12[1]/lbox[1])) * lbox[1];    

    r12[2] = r2[2] - r1[2];
    r12[2] -= static_cast<double>(Nint(r12[2]/lbox[2])) * lbox[2];
    
    return (SQ(r12[0]) + SQ(r12[1]) + SQ(r12[2]));
  }
  inline double distance3D(double *r12, const double *r1, const double *r2, const double *lbox,
			   const int *pbc_boundary){
    r12[0] = r2[0] - r1[0];
    r12[1] = r2[1] - r1[1];
    r12[2] = r2[2] - r1[2];
    
    if(pbc_boundary[0]) r12[0] -= static_cast<double>(Nint(r12[0]/lbox[0])) * lbox[0];
    if(pbc_boundary[1]) r12[1] -= static_cast<double>(Nint(r12[1]/lbox[1])) * lbox[1];
    if(pbc_boundary[2]) r12[2] -= static_cast<double>(Nint(r12[2]/lbox[2])) * lbox[2];
    return (SQ(r12[0])) + SQ(r12[1]) + SQ(r12[2]);
  }
  
  // place particle inside box: r_i in [0,L_i]
  inline void pbc(double &x, const double lbox){
    x = fmod(x + lbox, lbox);
    assert(x >= 0.0 && x < lbox);
  }
  inline void pbc(double *x, const double *lbox, const int &dim){
    for(int d = 0; d < dim; d++){
      x[d] = fmod(x[d] + lbox[d], lbox[d]);
      assert(x[d] >= 0.0 && x[d] < lbox[d]);
    }
  }
  inline void pbc(double *x, const double *lbox, const int *pbc_boundary, const int &dim){
    for(int d = 0; d < dim; d++){
      x[d] = (pbc_boundary[d] ? fmod(x[d] + lbox[d], lbox[d]) : x[d]);
      assert(!pbc_boundary[d] || (x[d] >= 0.0 && x[d] < lbox[d]));
    }
  }
  inline void pbc2D(double *x, const double *lbox){
    x[0] = fmod(x[0] + lbox[0], lbox[0]);
    x[1] = fmod(x[1] + lbox[1], lbox[1]);

    assert(x[0] >= 0.0 && x[0] < lbox[0]);
    assert(x[1] >= 0.0 && x[1] < lbox[1]);
  }
  inline void pbc2D(double *x, const double *lbox, const int *pbc_boundary){
    x[0] = (pbc_boundary[0] ? fmod(x[0] + lbox[0], lbox[0]) : x[0]);
    x[1] = (pbc_boundary[1] ? fmod(x[1] + lbox[1], lbox[1]) : x[1]);

    assert(!pbc_boundary[0] || (x[0] >= 0.0 && x[0] < lbox[0]));
    assert(!pbc_boundary[1] || (x[1] >= 0.0 && x[1] < lbox[1]));
  }
  inline void pbc3D(double *x, const double *lbox){
    x[0] = fmod(x[0] + lbox[0], lbox[0]);
    x[1] = fmod(x[1] + lbox[1], lbox[1]);
    x[2] = fmod(x[2] + lbox[2], lbox[2]);

    assert(x[0] >= 0.0 && x[0] < lbox[0]);
    assert(x[1] >= 0.0 && x[1] < lbox[1]);
    assert(x[2] >= 0.0 && x[2] < lbox[2]);
  }
  inline void pbc3D(double *x, const double *lbox, const int *pbc_boundary){
    x[0] = (pbc_boundary[0] ? fmod(x[0] + lbox[0], lbox[0]) : x[0]);
    x[1] = (pbc_boundary[1] ? fmod(x[1] + lbox[1], lbox[1]) : x[1]);
    x[2] = (pbc_boundary[2] ? fmod(x[2] + lbox[2], lbox[2]) : x[2]);

    assert(!pbc_boundary[0] || (x[0] >= 0.0 && x[0] < lbox[0]));
    assert(!pbc_boundary[1] || (x[1] >= 0.0 && x[1] < lbox[1]));
    assert(!pbc_boundary[2] || (x[2] >= 0.0 && x[2] < lbox[2]));
  }

  class LinkListOrtho{
  protected:
    static const int MAX_DIM  = 3;
    int Ncells, Ns[MAX_DIM];
    double iH[MAX_DIM];
  public:

    /*
      \brief Get optimum number of cells : cell size ~ rcut
     */
    virtual void get_opt_size(int *ns, const double &rcut, const double *lbox) const = 0;

    /*
      \brief Initialize link list for 2 or 3d simulations
     */
    virtual void init(int &num_cells, int &num_neighbors,    const int *ns, const double *lbox) = 0;

    /*
      Get unique cell coordinates for i-th cell
     */
    inline virtual void get_cellCoord(int *ixyz, const int &i) const = 0;
    
    /*
      Get list of cell coordinates and neighbor cell shifts
     */
    inline virtual void get_cells(int* cell, int* neighbors) const = 0;


    /*
      Get unique cell id given cell coordinates
     */
    inline virtual int get_cellID(const int *ixyz) const = 0;

    /*
      Get unique cell id given cell coordinates and neighbor shift
     */
    inline virtual int get_cellID(const int *ixyz, const int *dxyz) const = 0;

    /*
      Populate link list for given cell positions
     */
    virtual void populate_list(int *head, int *link,  double *r, const int &num) const = 0;

    /*
      Populate link list for given cell positions
     */
    virtual void populate_list(int *head, int *link,
			       double const* const* r, int const* stride, const int &num) const = 0;
  };
  
  class LinkListOrtho2D: public LinkListOrtho{
  protected:
    static const int DIM = 2;
    static const int NEIGHBORS = 5;
    static const int DNS[2*NEIGHBORS]; 
    
  public:
    LinkListOrtho2D(){
      Ncells = 0;
      Ns[0] = Ns[1] = Ns[2] = 1;
      iH[0] = iH[1] = iH[2] = 0.0;
    }

    void get_opt_size(int *ns, const double &rcut, const double *lbox) const {
      double delta = lbox[0];
      for(int d = 0; d < DIM; d++) delta = MIN(delta, lbox[d]);
      delta *= 1.0e-8;

      for(int d = 0; d < DIM; d++) ns[d] = static_cast<int>(lbox[d] / (rcut + delta));
    }

    void init(int &num_cells, int &num_neighbors,
	      const int *ns, const double *lbox){
      for(int d = 0; d < DIM; d++){
	run_assert(ns[d] > 0 && lbox[d] > 0.0, "Wrong link cell specs");	
	iH[d] = 1.0 / lbox[d];
	Ns[d] = ns[d];
      }
      
      {
	int too_small = 0;
	for(int d = 0; d < DIM; d++) if(Ns[d] < 3) too_small = 1;
	if(too_small) for(int d = 0; d < DIM; d++) Ns[d] = 1;
      }

      {
	Ncells = 1;
	for(int d = 0; d < DIM; d++) Ncells *= Ns[d];
      }
      run_assert(Ncells >= 1, "Ncells < 1");

      num_cells     = Ncells;
      num_neighbors = NEIGHBORS;
    }

    inline void get_cellCoord(int *ixyz, const int &i)const{
      ixyz[0] = static_cast<int>(i / Ns[1]);
      ixyz[1] = i - Ns[1]*ixyz[0];
    }
    
    inline void get_cells(int* cell, int* neighbors)const{
      int ii;
      if(cell != NULL){
	//set cell ids
	for(int i = 0; i < Ncells; i++){
	  ii = i * DIM;
	  this->get_cellCoord(&cell[ii], i);
	}
      }
      
      if(neighbors != NULL){
	for(int i = 0; i < NEIGHBORS; i++){
	  ii = i * DIM;
	  for(int d = 0; d < DIM; d++) neighbors[ii+d] = DNS[ii+d];
	}
      }
    }

    /*
      Get unique cell id given cell coordinates
     */
    inline int get_cellID(const int *ixyz)const{
      int ls[DIM] = {0, 0};
      ls[0] = (ixyz[0] + Ns[0]) % Ns[0];
      ls[1] = (ixyz[1] + Ns[1]) % Ns[1];      
      assert(ls[0] >= 0 && ls[0] < Ns[0]);
      assert(ls[1] >= 0 && ls[1] < Ns[1]);      
      return (Ns[1]*ls[0] + ls[1]);
    }

    /*
      Get unique cell id given cell coordinates and neighbor shift
     */
    inline int get_cellID(const int *ixyz, const int *dxyz)const{
      int ls[DIM] = {0, 0};
      ls[0] = (ixyz[0] + dxyz[0] + Ns[0]) % Ns[0];
      ls[1] = (ixyz[1] + dxyz[1] + Ns[1]) % Ns[1];      
      assert(ls[0] >= 0 && ls[0] < Ns[0]);
      assert(ls[1] >= 0 && ls[1] < Ns[1]);      
      return (Ns[1]*ls[0] + ls[1]);
    }

    /*
      Populate link list for given cell positions
     */
    void populate_list(int *head, int *link,
		       double *r, const int &num)const{
      for(int i = 0; i < Ncells; i++) head[i] = num;
      for(int i = 0; i < num; i++)    link[i] = num;
      
      int ls[DIM] = {0, 0};
      int ii, k;
      
      for(int i = 0; i < num; i++){
	
	ii = DIM*i;
	ls[0] = MIN(static_cast<int>((r[ii+0]*iH[0])*Ns[0]), Ns[0]-1);
	ls[1] = MIN(static_cast<int>((r[ii+1]*iH[1])*Ns[1]), Ns[1]-1);	
	k = (Ns[1]*ls[0] + ls[1]);
	
	if(head[k] == num){
	  head[k] = i;
	}else{
	  link[i] = head[k];
	  head[k] = i;
	}
	
      }//i
      
    }

    /*
      Populate link list for given cell positions
     */
    void populate_list(int *head, int *link,
		       double const* const* r, int const* stride, const int &num)const {
      for(int i = 0; i < Ncells; i++) head[i] = num;
      for(int i = 0; i < num; i++)    link[i] = num;
      
      int ls[DIM] = {0, 0};
      int k;
      
      for(int i = 0; i < num; i++){
	
	ls[0] = MIN(static_cast<int>((r[0][i*stride[0]]*iH[0])*Ns[0]), Ns[0]-1);
	ls[1] = MIN(static_cast<int>((r[1][i*stride[1]]*iH[1])*Ns[1]), Ns[1]-1);	
	k = (Ns[1]*ls[0] + ls[1]);
	
	if(head[k] == num){
	  head[k] = i;
	}else{
	  link[i] = head[k];
	  head[k] = i;
	}
      }// i
    }
  };

  /*const int LinkListOrtho2D::DNS[2*NEIGHBORS] = { 0,  0,	     
						  1,  0,
						  1,  1,
						  0,  1,
						  -1, 1 };*/
  
  /*const int LinkListOrtho::DNS3D[3*NEIGHBORS3D] = { 0,  0,  0,
						    1,  0,  0,
						    1,  1,  0,
						    0,  1,  0,
						    -1, 1,  0,
						    0,  0,  1,
						    1,  0,  1,
						    1,  1,  1,
						    0,  1,  1,
						    -1, 1,  1,
						    -1, 0,  1,
						    -1, -1, 1,
						    0,  -1, 1,
						    1,  -1, 1};*/

  
}



/*!
  Flag Waiving by Kevlin Henny
  C++ Workshop column in Application Development Advisor 6(3), April 2002
  http://www.two-sdg.demon.co.uk/curbralan/papers/FlagWaiving.pdf
  Usage:
  
  namespace STYLE{
    enum STYLE{OPT1, OPT2, OPT3, OPT4, ..., _SIZE_};
  }

  template<> struct enum_traits<STYLE::STYLE>: enum_traiter<STYLE::STYLE, STYLE::_SIZE>{};
  enum_set<STYLE::STYLE> SW_STYLE;  // this is the variable containing all the possible "STYLE" options
  
  if(SW_STYLE[STYLE::OPT1]) ...

*/
template<typename type>
struct enum_traits{
  static const bool        is_specialized = false;
  static const std::size_t count = 0;
};

template<typename type, type last_value>
struct enum_traiter{
  static const bool        is_specialized=true;
  static const std::size_t count = (last_value - type()) + 1;
};

template<typename enum_type,
	 typename traits = enum_traits<enum_type> >
class enum_set{
private:
  std::bitset<traits::count> bits;
  
public:
  enum_set(){
  }
  enum_set(enum_type setting){
    set(setting);
  }
  
  std::size_t count()const{
    return bits.count();
  }
  std::size_t size()const{
    return bits.size();
  }
  
  bool operator[](enum_type testing)const{
    assert(testing >= 0 && testing < traits::count);
    return bits[testing];
  }
  enum_set operator~()const{
    return enum_set(*this).flip();
  }
  bool any()const{
    return bits.any();
  }
  bool any(std::initializer_list<enum_type> testlist)const{
    bool res = false;
    for(auto &testing : testlist){
      assert(testing >= 0 && testing < traits::count);      
      res = (res || bits[testing]);
    }
    return res;
  }
  bool none()const{
    return bits.none();
  }
  bool none(std::initializer_list<enum_type> testlist)const{
    bool res = false;
    for(auto &testing : testlist){
      assert(testing >= 0 && testing < traits::count);            
      res = (res || bits[testing]);
    }      
    return !res;
  }
  bool all()const{
    return bits.all();
  }
  
  bool all(std::initializer_list<enum_type> testlist)const{
    bool res = true;
    for(auto &testing : testlist){
      assert(testing >= 0 && testing < traits::count);            
      res = (res && bits[testing]);
    }
    return res;
  }
  
  /*
    a  b  a AND b
    1  1   1
    1  0   0
    0  1   0
    0  0   0
  */
  bool any_and(std::initializer_list<enum_type> lista,
	       std::initializer_list<enum_type> listb){
    return (any(lista) && any(listb));
  }
  
  /*
    a  b  a XOR b
    1  1   0
    1  0   1
    0  1   1
    0  0   0
  */
  bool any_xor(std::initializer_list<enum_type> lista,
	       std::initializer_list<enum_type> listb){
    bool ta = any(lista);
    bool tb = any(listb);
    return (!(ta && tb) && (ta || tb));
  }
  
  /*
    a  b  a XNOR b (same as IF AND ONLY IF)
    1  1   1
    1  0   0
    0  1   0
    0  0   1
  */
  bool any_iff(std::initializer_list<enum_type> lista,
	       std::initializer_list<enum_type> listb){
    bool ta = any(lista);
    bool tb = any(listb);
    return ((ta && tb) || (!ta && !tb));
  }

  /*
    a  b  a OR b
    1  1   1
    1  0   1
    0  1   1
    0  0   0
  */
  bool any_or(std::initializer_list<enum_type> lista,
	      std::initializer_list<enum_type> listb){
    return (any(lista) || any(listb));
  }

  enum_set &operator&=(const enum_set &rhs){
    bits &= rhs.bits;
    return *this;
  }

  enum_set &operator|=(const enum_set &rhs){
    bits |= rhs.bits;
    return *this;
  }
  enum_set &operator ^=(const enum_set &rhs){
    bits ^= rhs.bits;
    return *this;
  }
  enum_set &set(){
    bits.set();
    return *this;
  }
  enum_set &set(enum_type setting, bool value=true){
    bits.set(setting, value);
    return *this;
  }
  enum_set &set(std::initializer_list<enum_type> setlist, bool value=true){
    for(auto &setting : setlist) bits.set(setting, value);
    return *this;
  }
  enum_set &reset(){
    bits.reset();
    return *this;
  }
  enum_set &reset(enum_type resetting){
    bits.reset(resetting);
    return *this;
  }
  enum_set &reset(std::initializer_list<enum_type> resetlist){
    for(auto &resetting : resetlist) bits.reset(resetting);
    return *this;
  }
  enum_set &flip(){
    bits.flip();
    return *this;
  }
  enum_set &flip(enum_type flipping){
    bits.flip(flipping);
    return *this;
  }
  enum_set &flip(std::initializer_list<enum_type> fliplist){
    for(auto &flipping : fliplist) bits.flip(flipping);
    return *this;
  }
};
#endif
