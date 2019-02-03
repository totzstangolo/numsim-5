
#include "typedef.hpp"
#include "geometry.hpp"
#include "parameter.hpp"
#include "grid.hpp"
#include "compute.hpp"
#include "iterator.hpp"
#include "solver.hpp"
#include "distri.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <unistd.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <memory>

using namespace std;

//extern unique_ptr<data> grid = make_unique<data>();

Data::Data(const index_t size){
	w = new real_t[9];
	index_x = new real_t[9];
	index_y = new real_t[9];
	f = new real_t* [18];
	for(int i = 0; i < 18; ++i){
		f[i] = new real_t[size];
	}
	boundary = new real_t[size];
	rho = new real_t[size];
	u = new real_t[size];
	v = new real_t[size];
}

Data::~Data(){
	if(rho!=nullptr){
		delete[] w;
		delete[] index_x;
		delete[] index_y;
		for(int i = 0; i < 18; ++i){
				delete[] f[i];
		}
		delete[] f;
		delete[] boundary;
		delete[] rho;
		delete[] u;
		delete[] v;
	}
}


Compute::Compute(const Geometry *geom, const Parameter *param){
	_geom = geom;
	_param = param;
	//debug init
	_max_dt = 1.0;

	// Initializing time, dtlimit and epslimit_
	_t = 0.0;
	multi_real_t h = _geom->Mesh();

	multi_real_t h2;
	h2[0] = h[0]*h[0];
	h2[1] = h[1]*h[1];
	_dtlimit = h2[0]*h2[1]*_param->Re()/(2*(h2[0]+h2[1]));

	if (_param->Pr()>0)
		_dtlimit = std::min<real_t>(_param->Re()*_param->Pr()/(2*(1/h2[0]+1/h2[1])),_dtlimit);

	_epslimit = _param->Eps() * _param->Eps() * _geom->Size()[0] * _geom->Size()[1];

	// creating grids with offset
	multi_real_t compute_offset;
	compute_offset[0] = -0.5 * h[0];
	compute_offset[1] = -0.5 * h[0];

	_f = new Distri(_geom,compute_offset);
	_f_new = new Distri(_geom,compute_offset);
	_f->Initialize(0.0);
	_f_new->Initialize(0.0);
	_f_eq = new Distri(_geom,compute_offset);
	_f_eq->Initialize(0.0);

	_u = new Grid(_geom, compute_offset);
	_u->Initialize(0);

	_v = new Grid(_geom, compute_offset);
	_v->Initialize(0);

	_p = new Grid(_geom, compute_offset);
	_T = new Grid(_geom, compute_offset);
	_tmp = new Grid(_geom, compute_offset);
	_p->Initialize(0);
	_T->Initialize(0);
	grid = new Data((_geom->Size()[0] + 2) * (_geom->Size()[1] + 2));
	Init();
	cout << _geom->Size()[0] + 2 << endl;


}

Compute::~Compute(){
	delete _u;
	delete _v;
	delete _p;
	// delete _T;
	// delete _rhs;
	// delete _tmp;
	// delete _solver;
	delete _f;
	delete _f_new;
}

void Compute::Init(){
	index_t n = (_geom->Size()[0] + 2) * (_geom->Size()[1] + 2);
	real_t tmp_w[9] = {4.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0, 1.0/9.0,
		  	  1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
	real_t tmp_index_x[9] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
	real_t tmp_index_y[9] = {0, 0, 0, 1, -1, 1, -1, -1, 1};
	for(index_t j = 0; j < 9; ++j){
		grid->w[j] = tmp_w[j];
		grid->index_x[j] = tmp_index_x[j];
		grid->index_y[j] = tmp_index_y[j];
	}

	for(index_t i = 0; i<n; ++i){
		 grid->boundary[i] = 0.0;
		 grid->rho[i] = 0.0;
		 grid->u[i] = 0.0;
		 grid->v[i] = 0.0;
		 for(index_t j = 0; j < 9; ++j) {
			 grid->f[j][i] = grid->w[j];
			 grid->f[j + 9][i] = grid->w[j];
		 }

	}
	  for(index_t i = 0; i < _geom->Size()[0] + 1; ++i){
		  grid->boundary[i] = 1;
		  grid->boundary[i * (_geom->Size()[0] + 2)] = 1;
		  grid->boundary[(i + 1) * (_geom->Size()[0] + 2) - 1] = 1;
		  //grid->boundary[i + n -(_geom->Size()[0] + 2)] = 1;
	  }

}

void Compute::TimeStep(bool printInfo){
	index_t n = (_geom->Size()[0] + 2) * (_geom->Size()[1] + 2);
	index_t m = _geom->Size()[0] + 2;
	//compute dt


	real_t dt = _param->Dt();
    if(dt == 0){
		dt = abs(std::min<real_t>(_geom->Mesh()[0]/_u->AbsMax(),_geom->Mesh()[1]/_v->AbsMax()));
		dt = std::min<real_t>(dt,_dtlimit);
		dt *= _param->Tau();
	}

	if(printInfo) {
		printf("dt: %f \n", dt);
	}
	/////////////////////////////////////////////////////////////////////////
	////////////// LATTICE BOLTZMANN IMPLEMENTATION (SECOND TRY) /////////////
	/////////////////////////////////////////////////////////////////////////

	// kinematic viscosity
    real_t nu    = _geom->GetInflowVelo()[0] *_geom->Size()[0] / _param->Re();
	real_t omega = 1.0 / (3.0*nu+1.0/2.0); //relaxation parameter

	///collision kernel for GPU
	// without iterator
	real_t cu,u2,v2, eq = 0.0;
	for(index_t i = 0; i < n; ++i){
		grid->rho[i] = 0;

		// computing rho, u and v
		for(index_t j = 0; j < 9; ++j){
			grid->rho[i] += grid->f[j][i];
		}
		grid->u[i] = (grid->f[1][i] - grid->f[2][i] + grid->f[5][i] - grid->f[6][i] + grid->f[7][i] - grid->f[8][i]) / grid->rho[i];
		grid->v[i] = (grid->f[3][i] - grid->f[4][i] + grid->f[5][i] - grid->f[6][i] - grid->f[7][i] + grid->f[8][i]) / grid->rho[i];
	}

	// inflow from north
	for(index_t i = n - m; i <n; ++i){
		/// rho
		grid->u[i] = _geom->GetInflowVelo()[0];
		grid->v[i] = _geom->GetInflowVelo()[1];
		grid->rho[i] = 1.0/(1.0 +_geom->GetInflowVelo()[1]) *
			(  grid->f[0][i] +  grid->f[1][i] + grid->f[2][i]
			 + 2.0 * (grid->f[3][i] + grid->f[5][i] + grid->f[8][i]));
		/// (Microscopic)
		grid->f[4][i] = grid->f[3][i] - 2/3 *  grid->rho[i] * _geom->GetInflowVelo()[1];
		grid->f[6][i] = grid->f[5][i] + 1/2*(grid->f[1][i] - grid->f[2][i])
					- 1/2*grid->rho[i]*_geom->GetInflowVelo()[0] - 1/6*grid->rho[i]*_geom->GetInflowVelo()[1];
		grid->f[7][i] = grid->f[8][i] + 1/2*(grid->f[2][i] - grid->f[1][i])
					+ 1/2*grid->rho[i]*_geom->GetInflowVelo()[0] - 1/6*grid->rho[i]*_geom->GetInflowVelo()[1];
	}


	for(index_t i = 0; i < n; ++i){
		//collision
		u2 = grid->u[i] * grid->u[i];
		v2 = grid->v[i] * grid->v[i];
		for(index_t j = 0; j < 9; ++j){
			cu = 3.0 * (grid->index_x[j] * grid->u[i] + grid->index_y[j] * grid->v[i]);
			eq = grid->rho[i] * grid->w[j] * (1 + cu + 0.5 * (cu * cu) - 1.5 * (u2 + v2));
			grid->f[9 + j][i] = grid->f[j][i] + omega * (eq - grid->f[j][i]);
		}
	}

	///streaming and bounce back
	int inv[9] = {0, 2, 1, 4, 3, 6, 5, 8, 7};
	int x_delta[9] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
	int y_delta[9] = {0, 0, 0, 1, -1, 1, -1, -1, 1};
	int k;
	//int debug;

	for(int h = 1; h < m - 1; ++h){
		for(int i = 1; i < m - 1; ++i){
			k = h * m + i;
			//cout <<k << " ";
			// Boundarys don't stream
			if(!grid->boundary[k]){
				grid->f[0][k] = grid->f[9][k];

				for(int j = 1; j < 9; ++j){
					//neighbor is boundary
					//debug = k - x_delta[j] - y_delta[j] * m;
					if(grid->boundary[k - x_delta[j] - y_delta[j] * m]){
						// bounceback
						grid->f[j][k] = grid->f[inv[j]+9][k];
					}else{
					//neighbor is no boundary (standard)
						grid->f[j][k] = grid->f[9 + j][k - x_delta[j] - y_delta[j] * m];
					}

				}

			}
		}

	}




	// iterator
	Iterator it(_geom);
	for(it.First();it.Valid();it.Next()){
		_f->rho(it) = _f->sum_vel(it);
		_p->Cell(it) = _f->rho(it);
		_u->Cell(it) = _f->sum_c_vel_x(it)/_f->rho(it);
		_v->Cell(it) = _f->sum_c_vel_y(it)/_f->rho(it);
	}

	BoundaryIterator lid(_geom);
	lid.SetBoundary(2);

	for(lid.First();lid.Valid();lid.Next()){
		/// Zhu, He (Macroscopic)
		_u->Cell(lid) = _geom->GetInflowVelo()[0];
		_v->Cell(lid) = _geom->GetInflowVelo()[1];
		_f->rho(lid) = 1.0/(1.0+_v->Cell(lid)) *
			(  _f->Cell(lid,0)+_f->Cell(lid,1)+_f->Cell(lid,3)
			 + 2*(_f->Cell(lid,2)+_f->Cell(lid,5)+_f->Cell(lid,6)));
		/// (Microscopic)
		_f->Cell(lid,4) = _f->Cell(lid,2) - 2/3*_f->rho(lid)*_v->Cell(lid);
		_f->Cell(lid,8) = _f->Cell(lid,6) + 1/2*(_f->Cell(lid,3)-_f->Cell(lid,1))
			+ 1/2*_f->rho(lid)*_u->Cell(lid) - 1/6*_f->rho(lid)*_v->Cell(lid);
		_f->Cell(lid,7) = _f->Cell(lid,5) + 1/2*(_f->Cell(lid,1)-_f->Cell(lid,3))
			- 1/2*_f->rho(lid)*_u->Cell(lid) - 1/6*_f->rho(lid)*_v->Cell(lid);
	}

	/// Collision
	for(it.First();it.Valid();it.Next()){
		u2 = _u->Cell(it)*_u->Cell(it);
		v2 = _v->Cell(it)*_v->Cell(it);
		for(int iCount=0;iCount<9;iCount++){
			cu = 3*(_f->ex(iCount)*_u->Cell(it)+_f->ey(iCount)*_v->Cell(it));
        	_f_eq->Cell(it,iCount) = _f->rho(it)*_f->t(iCount)
            * ( 1 + cu + 1/2*(cu*cu) - 3/2*(u2+v2));
        	_f_new->Cell(it,iCount) = _f->Cell(it,iCount)
				- omega * (_f->Cell(it,iCount) - _f_eq->Cell(it,iCount));
		}
	}

	// bounce back at walls
	real_t opp[9] = { 0, 3, 4, 1,  2, 7, 8,  5,  6};
	BoundaryIterator walls(_geom);
	walls.SetBoundary(0);
	for(walls.First();walls.Valid();walls.Next()){
		for (int iCount=0;iCount < 9;iCount++){
	        _f_new->Cell(walls,iCount) = _f->Cell(walls,opp[iCount]);
		}
	}
	walls.SetBoundary(1);
	for(walls.First();walls.Valid();walls.Next()){
		for (int iCount=0;iCount < 9;iCount++){
	        _f_new->Cell(walls,iCount) = _f->Cell(walls,opp[iCount]);
		}
	}
	walls.SetBoundary(3);
	for(walls.First();walls.Valid();walls.Next()){
		for (int iCount=0;iCount < 9;iCount++){
	        _f_new->Cell(walls,iCount) = _f->Cell(walls,opp[iCount]);
		}
	}

	// STREAMING
	InteriorIterator it2(_geom);
	for(it2.First();it2.Valid();it2.Next()){
		_f->Cell(it2,0) = _f_new->Cell(it2,0);
		_f->Cell(it2,1) = _f_new->Cell(it2.Right(),1);
		_f->Cell(it2,2) = _f_new->Cell(it2.Top(),2);
		_f->Cell(it2,3) = _f_new->Cell(it2.Left(),3);
		_f->Cell(it2,4) = _f_new->Cell(it2.Down(),4);
		_f->Cell(it2,5) = _f_new->Cell(it2.Top().Right(),5);
		_f->Cell(it2,6) = _f_new->Cell(it2.Top().Left(),6);
		_f->Cell(it2,7) = _f_new->Cell(it2.Down().Left(),7);
		_f->Cell(it2,8) = _f_new->Cell(it2.Down().Right(),8);
	}
	/////////////////////////////////////////////////////////////////////////
	////////////// END LATTICE BOLTZMANN IMPLEMENTATION /////////////////////
	/////////////////////////////////////////////////////////////////////////

	memcpy(&_u->Cell(0), &grid->u[0], sizeof(real_t) *n);
	memcpy(&_v->Cell(0), &grid->v[0], sizeof(real_t) *n);
	/*for(it.First();it.Valid();it.Next()){
		_u->Cell(it) = grid->u[it];
	}*/

	_t += dt;

	if(printInfo)
		// printf("time: %f \n itercount: %d \n max_dt: %f \n", _t, _iter_count, _max_dt);
		printf("time: %f\n", _t);

}

  /// Returns the simulated time in total
const real_t& Compute::GetTime() const{
    return _t;
}

  /// Returns the pointer to U
const Grid* Compute::GetU() const{
    return _u;
}
  /// Returns the pointer to V
const Grid* Compute::GetV() const{
    return _v;
}
  /// Returns the pointer to P
const Grid* Compute::GetP() const{
    return _p;
}

/// Returns the pointer to T
const Grid* Compute::GetT() const{
  return _T;
}

  /// Returns the pointer to RHS
const Grid* Compute::GetRHS() const{
    return _rhs;
}

/// Computes and returns the absolute velocity
const Grid* Compute::GetVelocity(){
	multi_real_t mid;
	Iterator iterator(_geom);
	for (iterator.First(); iterator.Valid(); iterator.Next()){
		mid[0] = _u->Cell(iterator) + _u->Cell(iterator.Left());
		mid[1] = _v->Cell(iterator) + _v->Cell(iterator.Down());
		_tmp->Cell(iterator) = sqrt((mid[0] * mid[0]) + (mid[1] * mid[1]));
	}
	return _tmp;
}

/// Computes and returns the vorticity
const Grid* Compute::GetVorticity(){
	// TODO real implementation needed
    Grid *grid = new Grid(_geom);
	return grid;
}

/// Computes and returns the stream line values
const Grid* Compute::GetStream(){
	// TODO real implementation needed
    Grid *grid = new Grid(_geom);
    return grid;
}

void Compute::NewVelocities(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_u->Cell(intIterator) = _F->Cell(intIterator) - (dt * _p->dx_r(intIterator));
		_v->Cell(intIterator) = _G->Cell(intIterator) - (dt * _p->dy_r(intIterator));
	}
}

/// Compute the temporary velocites F,G
void Compute::MomentumEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_F->Cell(intIterator) = _u->Cell(intIterator) +
			(
				 (_u->dxx(intIterator) + _u->dyy(intIterator))/ _param->Re()
				- _u->DC_udu_x(intIterator, _param->Alpha())
				- _u->DC_vdu_y(intIterator, _param->Alpha(), _v
			)) * dt;
		_G->Cell(intIterator) = _v->Cell(intIterator) +
			(
			 (_v->dyy(intIterator) + _v->dxx(intIterator))/ _param->Re()
				- _v->DC_vdv_y(intIterator, _param->Alpha())
				- _v->DC_udv_x(intIterator, _param->Alpha(), _u)
				) * dt;
	}

}

/// Compute the temporary velocites F,G
void Compute::ModMomentumEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	real_t tij = 0.0;
	real_t ti1j = 0.0;
	real_t tij1 = 0.0;
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		tij = _T->Cell(intIterator);
		ti1j = _T->Cell(intIterator.Right());
		tij1 = _T->Cell(intIterator.Top());
		_F->Cell(intIterator) = _F->Cell(intIterator) -
			dt*_param->Beta()*_param->Gx()*0.5*(tij+ti1j);
		_G->Cell(intIterator) = _G->Cell(intIterator) -
			dt*_param->Beta()*_param->Gy()*0.5*(tij+tij1);
	}
}

void Compute::TempEqu(const real_t &dt){
	InteriorIterator intIterator(_geom);
	real_t txx = 0.0;
	real_t tyy = 0.0;
	real_t utx = 0.0;
	real_t vty = 0.0;
	real_t pref = 1/(_param->Re()*_param->Pr());
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		txx = _T->dxx(intIterator);
		tyy = _T->dyy(intIterator);
		utx = _T->DC_duT_x(intIterator,_param->Alpha(),_u);
		vty = _T->DC_dvT_y(intIterator,_param->Alpha(),_v);
		_T->Cell(intIterator) = _T->Cell(intIterator) +
				(pref*(txx+tyy)-utx-vty)*dt;
		// std::cout << "Temperature in cell " << intIterator.Value() << ": "
		// 	<< _T->Cell(intIterator) << std::endl;
	}
	// exit(0);
}

/// Compute the RHS of the poisson equation
void Compute::RHS(const real_t &dt){
	InteriorIterator intIterator(_geom);
	for(intIterator.First(); intIterator.Valid(); intIterator.Next()){
		_rhs->Cell(intIterator) = (_F->dx_l(intIterator) + _G->dy_l(intIterator)) / dt;
	}
}
