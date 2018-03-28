#ifndef VICSEK_HXX
#define VICSEK_HXX

#include <cmath>
#include <exception>
#include <unistd.h>
#include "random.h"

#include "JsonParser.h"
#include "Common.hxx"
#include "Output.hxx"

typedef json_parser::json          json;

// link list parameters
MDutils::LinkListOrtho2D lnk;  // 2d orthorombic link list
int    lnk_cell_num;           // number of cells      
int    lnk_dcell_num;          // number of cell neighbors
nvec_i lnk_cell;               // cell ids
nvec_i lnk_dcell;              // cell neighbor shifts
vec_i  lnk_head;               // link list heads
vec_i  lnk_ihead;              // ordered heads
vec_i  lnk_list;               // link list

typedef void (*theta_frc0)(const int&i, const int& j, const double &r2, const double ni[DIM], const double nj[DIM],
			   double dni[DIM], double dnj[DIM]);
extern void theta_frc_sva(const int&, const int&, const double&, const double[], const double[],
			  double[], double[]);
extern void theta_frc_gca(const int&, const int&, const double&, const double[], const double[],
			  double[], double[]);
theta_frc0 theta_frc_ij;


void print_error_exit(const char* emsg){
  cerr << "# Error -> " << emsg << endl;
  exit(1);
}
template<typename T>
void print_error_exit(const char* emsg, const T &tag){
  cerr << "# Error -> " << emsg << " : " << tag << endl;
  exit(1);
}

// initalize simulation parameters
void read_input(const int &argc, char const* const* argv){
  using std::cerr;
  using std::endl;

  if(argc != 2) print_error_exit("Incorrect invocation", "\n# \tUsage : ./vicsek.x input_file.json");
  
  // Read input file and intialize parameter values
  string      string_flag;
  json        params;
  std::vector<int> dmy_ns(DIM, 1);  
  SW_NOISE.reset();
  SW_LINK.reset();
  SW_INIT.reset();
  SW_OUTPUT.reset();

  try{// parse input file
    json_parser::parse_file(argv[1], params);
    json_parser::dump("vicsek_params.out", params);
  }catch(exception &e){print_error_exit(e.what(), "parse input file");}
  
  try{// parse system parameters
    dir_name = json_parser::get<string>(params, "DIR");
    prj_name = json_parser::get<string>(params, "NAME");
    run_assert(dirCheckMake(dir_name.c_str()));

    { // Random number generator
      double seed;
      if((seed = json_parser::get<double>(params, "SEED")) > 0){
	RA_SEED(RND_SEED::GIVEN, seed);		
      }else{
	RA_SEED(RND_SEED::RANDOM, 0);	
      }
    }

    {// initial configuration
      string_flag = tolower(json_parser::get<string>(params, "INIT"));
      if(string_flag == "random"){
	SW_INIT.set(INIT::RANDOM);
      }else if(string_flag == "restart"){
	SW_INIT.set(INIT::RESTART);	  
      }else{
	print_error_exit("Unknown INIT option", string_flag);
      }
    }

    {// system parameters
      auto &loc = json_parser::get_child(params, "PARAMS");
      R0        = json_parser::get<double>(loc, "R0");
      DT        = json_parser::get<double>(loc, "DT");

      {	// adimensional buckingham PI parameters
	auto &ploc = json_parser::get_child(loc, "PI");
	PI_NOISE   = json_parser::get_vector<double>(ploc, "NOISE", 2);
	PI_VEL     = json_parser::get<double>(ploc, "VEL");
	PI_RHO     = json_parser::get<double>(ploc, "RHO");
	PI_BOX     = json_parser::get_vector<double>(ploc, "BOX", DIM);

	if(PI_NOISE[0] < 0.0 || PI_NOISE[0] > 1.0) print_error_exit("PI_NOISE[0] should lie in (0, 1)", PI_NOISE[0]);
	if(PI_NOISE[1] < 0.0 || PI_NOISE[1] > 1.0) print_error_exit("PI_NOISE[1] should lie in (0, 1)", PI_NOISE[1]);      
	if(PI_VEL < 0.0)  print_error_exit("PI_VEL cannot be negative", PI_VEL);
	if(PI_RHO <= 0.0) print_error_exit("PI_RHO should be positive", PI_RHO);
	if(PI_BOX[0] <= 0.0) print_error_exit("PI_BOX[0] should be positive", PI_BOX[0]);
	if(PI_BOX[1] <= 0.0) print_error_exit("PI_BOX[1] should be positive", PI_BOX[1]);
      }
	
      // noise              : ETA_ang, ETA_vec
      ETA_ANG    = PI_NOISE[0] / DT;
      ETA_VEC    = PI_NOISE[1] / DT;
      
      // interaction range
      RISQ    = SQ(R0);
      
      // box sizes          : Lx, Ly
      BOXL[0] = R0 / PI_BOX[0];
      BOXL[1] = R0 / PI_BOX[1];
      
      // densitiy
      RHO     = PI_RHO / (PI * RISQ);
      
      // number particles
      NUMP    = MAX(static_cast<int>(RHO * BOXL[0] * BOXL[1]), 1);
      
      // velocity
      U       = PI_VEL * (R0 / DT);
    }
  }catch(exception &e){print_error_exit(e.what(), "parse system parameters");}

  try{// parse output parameters
    auto &loc       = json_parser::get_child(params, "OUTPUT");
    NUM_FRAMES_SNAP = json_parser::get<int>(loc, "FRAMES_SNAP");
    NUM_SNAP        = json_parser::get<int>(loc, "NUMBER_SNAP");
    auto formats    = json_parser::get<std::vector<string>>(loc, "FORMAT");
    for(auto &fmt : formats){
      if(tolower(fmt) == "silo") SW_OUTPUT.set(OUTPUT::SILO);
      if(tolower(fmt) == "gsd")  SW_OUTPUT.set(OUTPUT::GSD);
      if(tolower(fmt) == "hdf5") SW_OUTPUT.set(OUTPUT::HDF5);
    }
    if(SW_OUTPUT.none({OUTPUT::SILO, OUTPUT::GSD, OUTPUT::HDF5})){
      print_error_exit("No valid output format specified : HDF5 | SILO | GSD");
    }
  }catch(exception &e){print_error_exit(e.what(), "parse output parameters");}
    
  try{// parse noise parameters
    auto &loc       = json_parser::get_child(params, "NOISE");
    {// ANGULAR noise
      if(json_parser::get<bool>(loc, "ANGULAR")){
	SW_NOISE.set(NOISE::ANGULAR);	
      }else{
	ETA_ANG = 0.0;
      }
    }
    {// VECTORIAL noise
      if(json_parser::get<bool>(loc["VECTORIAL"], "SELECT")){
	// Standard Vectorial Gregoire Chate Algorithm : Intrinsic random noise on i	
	SW_NOISE.set({NOISE::VECTORIAL});
	auto &vloc = json_parser::get_child(loc["VECTORIAL"], "ON");
	if(json_parser::get<bool>(vloc, "SELF_TERM")) SW_NOISE.set(NOISE::VECTORIAL_SELF); // SELF NOISE
      }else{
	ETA_VEC = 0.0;
      }
    }
  }catch(exception &e){print_error_exit(e.what(), "parse noise parameters");}
      
  try{// parse linked list parameters
    if(json_parser::get<bool>(params["LINKLIST"], "SELECT")){
      auto &loc = json_parser::get_child(params["LINKLIST"], "ON");
      if(json_parser::get<bool>(loc, "SORTED")) SW_LINK.set(LINK::SORTED);

      auto string_flag = json_parser::get<string>(loc, "SELECT", {"AUTO", "MANUAL"});
      if(string_flag == "AUTO"){
	SW_LINK.set({LINK::ON, LINK::AUTO});
	lnk.get_opt_size(dmy_ns.data(), sqrt(RISQ), BOXL);	
      }else if(string_flag == "MANUAL"){
	SW_LINK.set({LINK::ON, LINK::MANUAL});
	dmy_ns    = json_parser::get_vector<int>(loc, "NUMBER", DIM);
      }
    }else{
      SW_LINK.set(LINK::OFF);      
    }
    
    if(SW_LINK[LINK::ON]){
      lnk.init(lnk_cell_num, lnk_dcell_num, dmy_ns.data(), BOXL);
      if(lnk_cell_num == 1){
	cerr << "# Warning: Link Cell disabled! " << endl;
	SW_LINK.reset();
	SW_LINK.set(LINK::OFF);
      }
    }
    
  } catch(exception &e){print_error_exit(e.what(), "parse link list parameters");}
  
  run_assert(SW_NOISE.any(), "NOISE switch not set");
  run_assert(SW_LINK.any(), "LINK switch not set");
  run_assert(SW_INIT.any(), "INIT switch not set");
  run_assert(SW_OUTPUT.any(), "OUTPUT switch not set");

  run_assert(SW_NOISE.any_xor({NOISE::VECTORIAL, NOISE::VECTORIAL_SELF},
			      {NOISE::ANGULAR}),
	     "Noise VECTORIAL or ANGULAR flags should be set");
  run_assert(SW_NOISE.any_iff({NOISE::VECTORIAL_SELF},
			      {NOISE::VECTORIAL}),
	     "Noise VECTORIAL type should be set");
  
  
  run_assert(SW_LINK.any_xor({LINK::ON, LINK::AUTO, LINK::MANUAL, LINK::SORTED},
			     {LINK::OFF}),
	     "Link OFF or ON should be set");
  run_assert(SW_LINK.any_iff({LINK::ON}, {LINK::AUTO, LINK::MANUAL, LINK::SORTED}));
  
  run_assert(SW_INIT.any({INIT::RANDOM, INIT::RESTART}),
	     "Init RANDOM or RESTART should be set");  
  

  UDT = U*DT;
  MAX_NUM_STEPS = NUM_SNAP * NUM_FRAMES_SNAP;
  
  cerr << "# " << endl;
  cerr << "# Dir Name           = " << dir_name << endl;
  cerr << "# Prj Name           = " << prj_name << endl;
  cerr << "# " << endl;
  cerr << "# Output Formats     : " << endl;
  cerr << "# HDF5               = " << (SW_OUTPUT[OUTPUT::HDF5] ? "YES" : "NO") << endl;
  cerr << "# SILO               = " << (SW_OUTPUT[OUTPUT::SILO] ? "YES" : "NO") << endl;
  cerr << "# GSD                = " << (SW_OUTPUT[OUTPUT::GSD]  ? "YES" : "NO") << endl;
  cerr << "# " << endl;  
  cerr << "# Buckingham params  : " << endl;
  cerr << "# 1: PI_NOISE        = " << PI_NOISE[0] << " " << PI_NOISE[1] << endl;
  cerr << "# 2: PI_VEL          = " << PI_VEL << endl;
  cerr << "# 3: PI_RHO          = " << PI_RHO << endl;
  cerr << "# 4: PI_BOX          = " << PI_BOX[0]   << " " << PI_BOX[1] << endl;
  cerr << "# " << endl;
  cerr << "# System params      : " << endl;
  cerr << "# R0                 = " << R0 << endl;
  cerr << "# DT                 = " << DT << endl;
  cerr << "# Temp_ang/Temp_vec  = " << ETA_ANG << " " << ETA_VEC << endl;
  cerr << "# Speed              = " << U << endl;
  cerr << "# Density            = " << RHO << endl;  
  cerr << "# Num Particles      = " << NUMP << endl;
  cerr << "# Box Length         = " << BOXL[0] << " " << BOXL[1] << endl;
  cerr << "# R0_align           = " << sqrt(RISQ) << endl;
  cerr << "# Dump Steps         = " << NUM_FRAMES_SNAP << endl;
  cerr << "# Num Frames         = " << NUM_SNAP << endl;
  cerr << "# Max Steps          = " << MAX_NUM_STEPS << endl;
  cerr << "# " << endl;
  if(SW_LINK[LINK::ON]){
    cerr << "# Link Cell Ns       = ";
    for(int d = 0; d < DIM; d++) cerr << dmy_ns[d] << " ";
    cerr << endl;
    cerr << "# Cell Number        = " << lnk_cell_num << endl;
    cerr << "# Cell Neighbors     = " << lnk_dcell_num << endl;
  }
  cerr << "#"  << endl;

  {
    cerr << "# === OPTION SUMMARY === "  << endl;

    cerr << "# Noise Flags          : " << endl;
    cerr << "#   Angular  (SVA)     = " << (SW_NOISE[NOISE::ANGULAR] ? "ON" : "OFF") << endl;
    cerr << "#   Vectorial(GCA)     = " << (SW_NOISE[NOISE::VECTORIAL] ? "ON" : "OFF") << endl;
    if(SW_NOISE[NOISE::VECTORIAL]){
      cerr << "# \tself           = " << (SW_NOISE[NOISE::VECTORIAL_SELF] ? "ON" : "OFF") << endl;
    }
    cerr << "#" << endl;
    cerr << "# Link Flags           : " << endl;
    cerr << "#   Status             = " << (SW_LINK[LINK::OFF] ? "OFF" : "ON") << endl;
    if(SW_LINK[LINK::ON]){
      cerr << "# \tAuto           = " << (SW_LINK[LINK::AUTO] ? "ON" : "OFF")  << endl;
      cerr << "# \tManual         = " << (SW_LINK[LINK::MANUAL] ? "ON" : "OFF")<< endl;
      cerr << "# \tSorted         = " << (SW_LINK[LINK::SORTED] ? "ON" : "OFF")<< endl;
    }
    cerr << "#" << endl;
    cerr << "# Init Flags           : " << endl;
    cerr << "#   Random             = " << (SW_INIT[INIT::RANDOM] ? "ON" : "OFF") << endl;
    cerr << "#   Restart            = " << (SW_INIT[INIT::RESTART] ? "ON" : "OFF")<< endl;
    cerr << "#" << endl;
    cerr << "# === YRAMMUS NOITPO  === "  << endl;
    cerr << "#"  << endl;  
  }

}

inline void md_init(){
  init_threads();

  // allocate memory
  {
    Array::Allocate(theta, NUMP);
    Array::Allocate(xi, NUMP);
    
    Array::Allocate(sum_sin, NUMP);
    Array::Allocate(sum_cos, NUMP);
    Array::Allocate(sum_cnt, NUMP);    
    Array::Allocate(lop, NUMP);
    
    vxi.Allocate(NUMP, DIM);
    pos.Allocate(NUMP, DIM);
    vel.Allocate(NUMP, DIM);
    vlop.Allocate(NUMP, DIM);

    // zero data
    for(auto i = 0; i < NUMP; i++) theta[i]   = 0.0;
    for(auto i = 0; i < NUMP; i++) xi[i]      = 0.0;

    for(auto i = 0; i < NUMP; i++) sum_sin[i] = 0.0;
    for(auto i = 0; i < NUMP; i++) sum_cos[i] = 0.0;
    for(auto i = 0; i < NUMP; i++) sum_cnt[i] = 0;    
    for(auto i = 0; i < NUMP; i++) lop[i] = 0.0;
    vxi  = 0.0;
    pos  = 0.0;
    vel  = 0.0;
    vlop = 0.0;
  }
  
  // allocate link list memory
  if(lnk_cell_num > 1){
    lnk_dcell.Allocate(lnk_dcell_num, DIM);
    lnk.get_cells(NULL, &lnk_dcell[0][0]);    

    Array::Allocate(lnk_head, lnk_cell_num);
    Array::Allocate(lnk_list, NUMP);
    if(SW_LINK[LINK::SORTED]) Array::Allocate(lnk_ihead, lnk_cell_num);
  }

  // force routines
  {
    if(SW_NOISE[NOISE::VECTORIAL]){
      theta_frc_ij = theta_frc_gca;
    }else{
      theta_frc_ij = theta_frc_sva;
    }
  }
}

inline void configuration_init(){
  // Initial configuration
  if(SW_INIT[INIT::RANDOM]){
    cerr << "# Initial configuration is RANDOM" << endl;
    for(int i = 0; i < NUMP; i++){
      theta[i] = RA_CO(-PI, PI);
      
      vec ri = pos[i];
      ri[0]  = RA_CO(BOXL[0]);
      ri[1]  = RA_CO(BOXL[1]);	
    }
  }else if(SW_INIT[INIT::RESTART]){
    string restart_file = prj_name + ".restart.h5";
    cerr << "# Initial configuration is RESTART: " << restart_file <<endl;      
    h5out.open_file(restart_file.c_str(), "r");
    h5out.read_data("pos", h5vector_space, &pos[0][0]);
    h5out.read_data("theta", h5scalar_space, &theta[0]);
    h5out.close_file();
  }
  
  for(auto i = 0; i < NUMP; i++){
    vec vi = vel[i];
    vi[0]  = cos(theta[i]);
    vi[1]  = sin(theta[i]);
  }
}

inline void write_configuration(){
  static int frame_id = 0;
  if(SW_OUTPUT[OUTPUT::HDF5]) hdf5_output(frame_id);
  if(SW_OUTPUT[OUTPUT::SILO]) silo_output(frame_id);
  if(SW_OUTPUT[OUTPUT::GSD])  gsd_output(frame_id);
  frame_id++;
}

inline void output_init(){
  // Dataspaces for IO
  { 
    SimIO::io_size rank, dims[DIM];    
    // dataspace for scalar data
    rank = 1;
    dims[0] = NUMP;
    h5scalar_space = SimIO::Space::create(rank, dims);
    
    // dataspace for vector data
    rank = 2;
    dims[0] = NUMP;
    dims[1] = DIM;
    h5vector_space = SimIO::Space::create(rank, dims);
  }
  if(SW_OUTPUT[OUTPUT::HDF5]) hdf5_init();
  if(SW_OUTPUT[OUTPUT::SILO]) silo_init();
  if(SW_OUTPUT[OUTPUT::GSD])  gsd_init();
  if(SW_OUTPUT.any({OUTPUT::SILO, OUTPUT::GSD})){
    float_buffer.Allocate(NUMP*3);
    uint_buffer.Allocate(NUMP);
  }
}

inline void output_end(){

  if(SW_OUTPUT[OUTPUT::HDF5]) hdf5_end();
  if(SW_OUTPUT[OUTPUT::SILO]) silo_end();
  if(SW_OUTPUT[OUTPUT::GSD])  gsd_end();
  
  {// write restart file
    h5out.create_file((prj_name+".restart.h5").c_str());
    h5out.write_data("theta", h5scalar_space, &theta[0], static_cast<double*>(&theta[0]));
    h5out.write_data("pos",   h5vector_space, &pos[0][0],static_cast<double*>(&pos[0][0]));
    h5out.close_file();
  }
  {// close dataspaces
    SimIO::Space::close(h5scalar_space);
    SimIO::Space::close(h5vector_space);
  }
}

inline void initialize(int argc, char** argv){
  read_input(argc, argv);
  md_init();
  output_init();
  configuration_init();
}
inline void finalize(){
  output_end();
}

#endif
