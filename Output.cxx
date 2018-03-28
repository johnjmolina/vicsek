#include "Output.hxx"

string dir_name;
string prj_name;

SimIO::SimIO             h5out;   // hdf5 output writer
SimIO::io_space h5scalar_space;  // dataspace for 1d data
SimIO::io_space h5vector_space;  // dataspace for 2d data

gsd_handle      gsdout;
Array::array1<float>         float_buffer;
Array::array1<unsigned int>  uint_buffer;

void hdf5_init(){
  h5out.create_file((dir_name + "/" + prj_name + ".h5").c_str());
  {
    static const char sys_name[8] = "/sys";
    static const char trj_name[8] = "/trj";
    static const double eta[2] = {ETA_ANG, ETA_VEC};
    h5out.create_group(sys_name);      
    h5out.write_attr("dt", DT);
    h5out.write_attr("eta", eta, 2);
    h5out.write_attr("rho", RHO);
    h5out.write_attr("R0", R0);
    h5out.write_attr("u", U);
    h5out.write_attr("num", NUMP);  
    h5out.write_attr("box", BOXL, 2);
    h5out.write_attr("frames", NUM_SNAP);
    h5out.write_attr("angular_noise", SW_NOISE[NOISE::ANGULAR] ? "on" : "off");
    
    if(!SW_NOISE[NOISE::VECTORIAL]){
      h5out.write_attr("vectorial_noise", "off");      	
    }
    
    {
      SimIO::io_size rank, dims[DIM];      
      rank = 2;
      dims[0] = DIM;
      dims[1] = DIM;
      auto dmy_space = SimIO::Space::create(rank, dims);
      double xx[DIM*DIM] = {0.0, 0.0, BOXL[0], BOXL[0]};
      double yy[DIM*DIM] = {0.0, BOXL[1], 0.0, BOXL[1]};
      h5out.write_data("box_x", dmy_space, xx);
      h5out.write_data("box_y", dmy_space, yy);
      SimIO::Space::close(dmy_space);      
    }
    h5out.close_group(); // close sys group
    h5out.create_group(trj_name); // open trj group
  }  
}

void hdf5_output(int frame_id){
  char frame_name[BUFFER_SIZE];  
  sprintf(frame_name, "./t_%d", frame_id);
  h5out.create_group(frame_name);
  h5out.write_attr("ts", ts);
  h5out.write_data("num", h5scalar_space, &sum_cnt[0]);
  h5out.write_data("phi", h5scalar_space, &lop[0]);
  h5out.write_data("pol", h5vector_space, &vlop[0][0]);  
  h5out.write_data("pos", h5vector_space, &pos[0][0]);
  h5out.write_data("vel", h5vector_space, &vel[0][0]);
  h5out.close_group();
}

void hdf5_end(){
  //close trajectory file
  h5out.close_file();  
  
  {
    char buffer[BUFFER_SIZE];
    snprintf(buffer, BUFFER_SIZE,
	     "./vicsek_xdmf.py -d %s -f %s.h5 -t %d -n %d",
	     dir_name.c_str(), prj_name.c_str(), NUM_SNAP + 1, NUMP
	     );
    cerr << "# To generate xdmf file run:" << endl;  
    cerr << "# " << buffer << endl;  
    //    system(buffer);
    cerr << "# The end!" << endl;
  }
}

void silo_init(){
}
void silo_output(int frame_id){
  static string silo_filename = ("./" + dir_name + "/" + prj_name);
  static const char* vnames[] = {"v_x", "v_y"};  
  DBfile *dbfile = NULL;    
  char buffer[BUFFER_SIZE];
  double time;
  snprintf(buffer, BUFFER_SIZE, "%s_%06d.silo", silo_filename.c_str(), frame_id);
  dbfile = DBCreate(buffer, DB_CLOBBER, DB_LOCAL, "vicsek", DB_HDF5);
  DBoptlist *optlist = DBMakeOptlist(2);

  Array::array2<float> tvec(DIM, NUMP, float_buffer());  
  {
    time = static_cast<double>(ts);
    DBAddOption(optlist, DBOPT_DTIME, &time);
    DBAddOption(optlist, DBOPT_CYCLE, &frame_id);
  }
  {
    double rx[] = {0.0, BOXL[0]};
    double ry[] = {0.0, BOXL[1]};
    double* coords[] = {rx, ry};
    int dims[]  = {2,2};
    DBPutQuadmesh(dbfile, "box", NULL, coords, dims, DIM, DB_DOUBLE, DB_COLLINEAR, optlist);
  }
  {
    int dims[] = {NUMP};

#pragma omp parallel for
    for(int i = 0; i < NUMP; i++) tvec(0,i) = pos(i,0);
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++) tvec(1,i) = pos(i,1);
    float * coords[DIM] = {tvec[0], tvec[1]};
    DBPutPointmesh(dbfile, "particles", DIM, coords, NUMP, DB_FLOAT, NULL);
    DBPutQuadvar1(dbfile, "num",  "particles", &sum_cnt[0], dims, 1, NULL, 0, DB_INT, DB_NODECENT, NULL);
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++) tvec(0,i) = lop[i];
    DBPutQuadvar1(dbfile, "phi", "particles", tvec[0], dims, 1, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
  }
  {
    int dims[] = {NUMP};
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++) tvec(0,i) = U*vel(i,0);
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++) tvec(1,i) = U*vel(i,1);
    float* v_ptr[DIM] = {tvec[0], tvec[1]};
    DBPutQuadvar(dbfile, "vel", "particles", DIM, vnames, v_ptr, dims, 1, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
  }
  DBFreeOptlist(optlist);
  DBClose(dbfile);  
}
void  silo_end(){
}

template<typename T, typename Q>
inline void gsdCopyVector(const Array::array2<T> orig2d, const T *shift, const T &scale, const T &zval, Q* ptr){
  Array::array2<Q> dmy3d(NUMP, 3, ptr);
  for(auto i = 0; i < NUMP; i++){
    dmy3d(i, 0) = (orig2d(i,0) - shift[0])*scale;
    dmy3d(i, 1) = (orig2d(i,1) - shift[1])*scale;
    dmy3d(i, 2) = zval;
  }
}
void gsd_init(){
  const char* filename = (dir_name + "/" + prj_name + ".gsd").c_str();
  run_assert(gsd_create(filename, "Vicsek", "hoomd", gsd_make_version(1,1)) == 0, \
	     "gsd create");
  run_assert(gsd_open(&gsdout, filename, GSD_OPEN_APPEND) == 0, \
	     "gsd open");
}
void gsd_output(int frame_id){
  {
    uint64_t _ts = ts;
    uint8_t  dim= DIM;
    float    box[6] = {static_cast<float>(BOXL[0]),
		       static_cast<float>(BOXL[1]),
		       1.0,
		       0.0, 0.0, 0.0};    
    
    run_assert(gsd_write_chunk(&gsdout, "configuration/step", GSD_TYPE_UINT64, 1, 1, 0, &_ts) == 0,
	       "gsd write configuration/step");
    run_assert(gsd_write_chunk(&gsdout, "configuration/dimensions", GSD_TYPE_UINT8, 1, 1, 0, &dim)==0,
	       "gsd write configuration/dimensions");
    run_assert(gsd_write_chunk(&gsdout, "configuration/box", GSD_TYPE_FLOAT, 6, 1, 0, &box) == 0,
	       "gsd write configuration/box");
  }

  { 
    {// particle data
      uint32_t n = NUMP;
      run_assert(gsd_write_chunk(&gsdout, "particles/N", GSD_TYPE_UINT32, 1, 1, 0, &n) == 0,
		 "gsd write particles/N");
    }

    {// constant data
      if(ts == 0){
#pragma omp parallel for
	for(auto i = 0; i < NUMP; i++) float_buffer(i) = 0.1*R0;
	run_assert(gsd_write_chunk(&gsdout, "particles/diameter", GSD_TYPE_FLOAT, NUMP, 1, 0, float_buffer) == 0,
		   "gsd write particles/diameter");
      }
    }
    
    {// positions
      double shift[DIM] = {BOXL[0]/2.0, BOXL[1]/2.0};    
      gsdCopyVector(pos, shift, 1.0, 0.01*R0, float_buffer());
      run_assert(gsd_write_chunk(&gsdout, "particles/position", GSD_TYPE_FLOAT, NUMP, 3, 0, float_buffer) == 0,
		 "gsd write particles/position");
    }
    
    {// velocities
      double shift[DIM] = {0.0, 0.0};
      gsdCopyVector(vel, shift, 2*R0, 0.0, float_buffer());
      run_assert(gsd_write_chunk(&gsdout, "particles/velocity", GSD_TYPE_FLOAT, NUMP, 3, 0, float_buffer) == 0,
		 "gsd write particles/velocity");
    }
  }
  {
    run_assert(gsd_end_frame(&gsdout)== 0, "gsd write frame");
  }
}
void gsd_end(){
  gsd_close(&gsdout);
}
