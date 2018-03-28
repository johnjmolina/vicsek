#include "Vicsek.hxx"
inline void theta_frc_sva(const int &i, const int &j, const double &r2, const double ni[DIM], const double nj[DIM],
			  double dni[DIM], double dnj[DIM]){
  dni[0] = nj[0];
  dni[1] = nj[1];
  
  dnj[0] = ni[0];
  dnj[1] = ni[1];
}
inline void theta_frc_gca(const int &i, const int &j, const double &r2, const double ni[DIM], const double nj[DIM],
			  double dni[DIM], double dnj[DIM]){
  if( i != j){
    vec csxi = vxi[i];      
    dni[0] = nj[0] + ETA_VEC*csxi[0];
    dni[1] = nj[1] + ETA_VEC*csxi[1];
    
    vec csxj = vxi[j];      
    dnj[0]   = ni[0] + ETA_VEC*csxj[0];
    dnj[1]   = ni[1] + ETA_VEC*csxj[1];
  }else{ // i = j 
    vec csxi      = vxi[i];
    dni[0] = dnj[0] = (ni[0] + ETA_VEC*csxi[0]);
    dni[1] = dnj[1] = (ni[1] + ETA_VEC*csxi[1]);
  }
}

/*
  compute pair interactions : average sin / cos 
*/
inline void md_force_solver(){

#pragma omp parallel for
  for(int i = 0; i < NUMP; i++) sum_sin[i] = 0.0;
#pragma omp parallel for  
  for(int i = 0; i < NUMP; i++) sum_cos[i] = 0.0;
#pragma omp parallel for  
  for(int i = 0; i < NUMP; i++) sum_cnt[i] = 0;
#pragma omp parallel for
  for(int i = 0; i < NUMP; i++) vlop[i][0] = vlop[i][1] = 0.0;

  //precompute vectorial noise
  if(SW_NOISE[NOISE::VECTORIAL]){
    for(int i = 0; i < NUMP; i++){
      vec csxi = vxi[i];
      xi[i] = RA_CO(-PI, PI);
      csxi[0] = cos(xi[i]);
      csxi[1] = sin(xi[i]);
    }
  }

  //compute link lists
  if(SW_LINK[LINK::ON]){
    lnk.populate_list(lnk_head, lnk_list, pos, NUMP);
    if(SW_LINK[LINK::SORTED]){
      for(int i = 0; i < lnk_cell_num; i++) lnk_ihead[i] = i;      	
      std::sort(&lnk_ihead[0],
		&lnk_ihead[0] + lnk_cell_num,
		[](const int &a, const int &b){return lnk_head[a] < lnk_head[b];});
    }

    int max_iic  = lnk_cell_num;
    if(SW_LINK[LINK::SORTED]) for(max_iic=0; max_iic < lnk_cell_num && lnk_head[lnk_ihead[max_iic]] < NUMP; max_iic++);
    
#pragma omp parallel for
    for(int iic=0; iic < max_iic; iic++){ // over cell ids
      int ic = (SW_LINK[LINK::SORTED] ? lnk_ihead[iic] : iic);
      int icell[DIM];      
      double r2, rij[DIM], dni[DIM], dnj[DIM];
      
      for(int inc=0; inc < lnk_dcell_num; inc++){ // over neighboring cells

	lnk.get_cellCoord(icell, ic);
	int jc = lnk.get_cellID(icell, lnk_dcell[inc]); // neighbor cell id
	//int jc = lnk.get_cellID(lnk_cell[ic], lnk_dcell[inc]);
	
	int i = lnk_head[ic];
	
	while(i < NUMP){ // over i particles in center cell
	  vec vi = vel[i];	  
	  int j = lnk_head[jc];
	  while(j < NUMP){ // over j particles in neighbor cell
	    // valid pairs withough double counting
	    if((ic != jc || j >= i) && ((r2 = distance2D(rij, pos[i], pos[j], BOXL)) <= RISQ)){
	      vec vj = vel[j];
	      theta_frc_ij(i, j, r2, vi, vj, dni, dnj);

#pragma omp atomic
	      sum_cos[i] += dni[0];
#pragma omp atomic
	      sum_sin[i] += dni[1];
#pragma omp atomic
	      vlop[i][0] += vj[0];
#pragma omp atomic
	      vlop[i][1] += vj[1];
#pragma omp atomic
	      sum_cnt[i] += 1;
	      
	      if(i != j){
#pragma omp atomic
		sum_cos[j] += dnj[0];
#pragma omp atomic
		sum_sin[j] += dnj[1];
#pragma omp atomic
		vlop[j][0] += vi[0];
#pragma omp atomic
		vlop[j][1] += vi[1];
#pragma omp atomic
		sum_cnt[j] += 1;
	      } // i != j
	      
	    }// valid pairs
	    j = lnk_list[j];
	    
	  }// j < NUMP
	  i = lnk_list[i];
	  
	}// i < NUMP
	
      }// inc
      
    }// iic
    
  }else{
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++){
      int cnt_align, cnt_noise;
      double r2, rij[DIM], dni[DIM], dnj[DIM];
      vec vi = vel[i];
      
      for(int j = i; j < NUMP; j++){
	if( (r2 = distance2D(rij, pos[i], pos[j], BOXL)) <= RISQ){
	  vec vj = vel[j];
	  theta_frc_ij(i, j, r2, vi, vj, dni, dnj);

#pragma omp atomic
	  sum_cos[i] += dni[0];
#pragma omp atomic
	  sum_sin[i] += dni[1];
#pragma omp atomic
	  vlop[i][0] += vj[0];
#pragma omp atomic
	  vlop[i][1] += vj[1];
#pragma omp atomic
	  sum_cnt[i] += 1;

	  if(i != j){
#pragma omp atomic
	    sum_cos[j] += dnj[0];
#pragma omp atomic
	    sum_sin[j] += dnj[1];
#pragma omp atomic
	    vlop[j][0] += vi[0];
#pragma omp atomic
	    vlop[j][1] += vi[1];
#pragma omp atomic
	    sum_cnt[j] += 1;
	  }
	}
      }//j
    }//i
  }
#pragma omp parallel for
  for(int i = 0; i < NUMP; i++) lop[i] = sqrt(SQ(vlop[i][0]) + SQ(vlop[i][1])) / static_cast<double>(sum_cnt[i]);
}
inline void md_theta_solver(){
  /*
    Compute new orientation angle
  */
  if(SW_NOISE[NOISE::ANGULAR]){ // YES ANGULAR NOISE
    for(int i = 0; i < NUMP; i++) xi[i] = RA_CO(-ETA_ANG*PI, ETA_ANG*PI);        
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++){
      theta[i] = atan2(sum_sin[i], sum_cos[i]) + xi[i];      
    }
  }else{ // NO ANGULAR NOISE
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++){
      theta[i] = atan2(sum_sin[i], sum_cos[i]);
    }
  }
}
inline void md_pos_solver(){
  /*
    Update positions & velocities
  */
  {
#pragma omp parallel for
    for(int i = 0; i < NUMP; i++){
      vec ri = pos[i];
      vec vi = vel[i];

      vi[0]  =  cos(theta[i]);
      vi[1]  =  sin(theta[i]);
      
      ri[0] +=  vi[0] * UDT;//x
      ri[1] +=  vi[1] * UDT;//y
      pbc2D(ri, BOXL);
    }
  }
}

inline void md_step(const bool &dump){
  /*
    Compute pair interactions and update particle orientations & positions
   */
  md_force_solver();
  if(dump) write_configuration();
  
  md_theta_solver();
  md_pos_solver();
}

int main(int argc, char** argv){
  initialize(argc, argv);
  WallTimer timer;
  timer.start();
  int CHECK_POINT = MAX(MAX_NUM_STEPS / 10, 1);
  cerr << "# 0  10  20  30  40  50  60  70  80  90  100%\n# *";
  for(ts = 0; ts < MAX_NUM_STEPS; ts++){
    md_step((ts % NUM_FRAMES_SNAP) == 0);        
    if(ts % CHECK_POINT == 0 && ts > 0) cerr << "****";
  }
  if(ts % NUM_FRAMES_SNAP == 0){
    md_force_solver();
    write_configuration();
  }
  cerr << "****" << endl;
  double wall_time = timer.stop();
  cerr << "# Execution Finished" << endl;
  cerr << "# \tSeconds : " << wall_time << endl;
  cerr << "# \tMinutes : " << wall_time / 60.0 << endl;
  cerr << "# \tHours   : " << wall_time / 3600.0 << endl;
  cerr << "# \tDays    : " << wall_time / 86400.0 << endl;

  finalize();
  return(0);
}
