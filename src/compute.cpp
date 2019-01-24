
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
using namespace std;

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
	_f->Initialize(1.0);
	_f_new->Initialize(1.0);

	_u = new Grid(_geom, compute_offset);
	_u->Initialize(0);
	BoundaryIterator itBound(_geom);
	itBound.SetBoundary(2);
	for(itBound.First();itBound.Valid();itBound.Next()){
		// _f_new->Cell(itBound.Down(),8) += _geom->GetInflowVelo()[0];
		// _f_new->Cell(itBound.Down(),7) -= _geom->GetInflowVelo()[0];
		_u->Cell(itBound) = _geom->GetInflowVelo()[0];
	}

	_v = new Grid(_geom, compute_offset);
	_v->Initialize(0);

	_p = new Grid(_geom, compute_offset);
	_T = new Grid(_geom, compute_offset);
	_tmp = new Grid(_geom, compute_offset);
	_p->Initialize(0);
	_T->Initialize(0);
}

Compute::~Compute(){
	delete _u;
	delete _F;
	delete _v;
	delete _G;
	delete _p;
	delete _T;
	delete _rhs;
	delete _tmp;
	delete _solver;
	delete _f;
	delete _f_new;
}

void Compute::TimeStep(bool printInfo){
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
	////////////// LATTICE BOLTZMANN IMPLEMENTATION (FIRST TRY) /////////////
	/////////////////////////////////////////////////////////////////////////
	/// COLLISION

	InteriorIterator it2(_geom);
	real_t r=0;
	real_t fe[9];
	real_t vis = 0.01;
	real_t omega = 1.0/(0.5 + 3.0 * vis);
	real_t eomega = 1.0-omega;
	real_t tau = 0.0001;
	real_t u2 = 0;
	real_t uTot = 0.0;
	real_t uTot2 = 0.0;
	real_t e[18] = {0.0,0.0,  1.0,0.0,  0.0,1.0,
		 		   -1.0,0.0,  0.0,-1.0, 1.0,1.0,
	   			   -1.0,1.0, -1.0,-1.0, 1.0,-1.0};
	real_t ti[9] = {4.0/9.0,  1.0/9.0,  1.0/9.0,   1.0/9.0,   1.0/9.0,
							  1.0/36.0, 1.0/36.0,  1.0/36.0,  1.0/36.0};
	for(it2.First();it2.Valid();it2.Next()){
		r = 0;
		_u->Cell(it2) = _f_new->Cell(it2,1)*e[2] + _f_new->Cell(it2,3)*e[6] +
						_f_new->Cell(it2,5)*e[10] + _f_new->Cell(it2,6)*e[12]
						+ _f_new->Cell(it2,7)*e[14] + _f_new->Cell(it2,8)*e[16];
		_v->Cell(it2) = _f_new->Cell(it2,2)*e[5] + _f_new->Cell(it2,4)*e[9] +
						_f_new->Cell(it2,5)*e[11] + _f_new->Cell(it2,6)*e[13]
						+ _f_new->Cell(it2,7)*e[15] + _f_new->Cell(it2,8)*e[17];
		for (index_t i=0;i<9;i++)
			r += _f_new->Cell(it2,i);
		_u->Cell(it2) /= r;
		if(abs(_t)>1.115){exit(1);}
		_v->Cell(it2) /= r;
		_p->Cell(it2) = r;

		// e[18] = {0.0,0.0,  1.0,0.0,  0.0,1.0,
		// 		   -1.0,0.0,  0.0,-1.0, 1.0,1.0,
		// 		   -1.0,1.0, -1.0,-1.0, 1.0,-1.0};
		for (int iCount=0;iCount<18;iCount+=2){
			real_t uTot = e[iCount]*_u->Cell(it2)+e[iCount+1]*_v->Cell(it2);
			real_t uTot2 = uTot*uTot;
			fe[iCount/2] = ti[iCount/2]*r*(1.0 + 3*uTot + 4.5*uTot2
				- 1.5*(_u->Cell(it2)*_u->Cell(it2) + _v->Cell(it2)*_v->Cell(it2)));

		}
		printf ("F_0^eq: %f, Density: %f\n",fe[1],r);
		for (int iCount=0;iCount<9;iCount++){
			_f->Cell(it2,iCount) = omega*dt*fe[iCount] + eomega*dt*_f_new->Cell(it2,iCount);
		}
	}

	/// STREAMING
	BoundaryIterator itBound(_geom);
	// itBound.SetBoundary(2);
	// for(itBound.First();itBound.Valid();itBound.Next()){
	// 	_f_new->Cell(itBound,1) += _geom->GetInflowVelo()[0];
	// 	// _f->Cell(itBound.Down(),8) += _geom->GetInflowVelo()[0];
	// }
	InteriorIterator it(_geom);
	for(it.First();it.Valid();it.Next()){
		_f_new->Cell(it.Right(),1) = _f->Cell(it,1);
		_f_new->Cell(it.Top(),2) = _f->Cell(it,2);
		_f_new->Cell(it.Left(),3) = _f->Cell(it,3);
		_f_new->Cell(it.Down(),4) = _f->Cell(it,4);
		_f_new->Cell(it.Right().Top(),5) = _f->Cell(it,5);
		_f_new->Cell(it.Left().Top(),6) = _f->Cell(it,6);
		_f_new->Cell(it.Left().Down(),7) = _f->Cell(it,7);
		_f_new->Cell(it.Right().Down(),8) = _f->Cell(it,8);
		_f_new->Cell(it,0) = _f->Cell(it,0);
	}
	// BoundaryIterator itBound(_geom);
	itBound.SetBoundary(0);
	for(itBound.First();itBound.Valid();itBound.Next()){
		_f_new->Cell(itBound,2) = _f_new->Cell(itBound,4);
		_f_new->Cell(itBound,5) = _f_new->Cell(itBound,7);
		_f_new->Cell(itBound,6) = _f_new->Cell(itBound,8);
	}
	itBound.SetBoundary(1);
	for(itBound.First();itBound.Valid();itBound.Next()){
		_f_new->Cell(itBound,1) = _f_new->Cell(itBound,3);
		_f_new->Cell(itBound,5) = _f_new->Cell(itBound,7);
		_f_new->Cell(itBound,8) = _f_new->Cell(itBound,6);
	}
	itBound.SetBoundary(3);
	for(itBound.First();itBound.Valid();itBound.Next()){
		_f_new->Cell(itBound,3) = _f_new->Cell(itBound,1);
		_f_new->Cell(itBound,6) = _f_new->Cell(itBound,8);
		_f_new->Cell(itBound,7) = _f_new->Cell(itBound,5);
	}
	// itBound.SetBoundary(2);
	// for(itBound.First();itBound.Valid();itBound.Next()){
	// 	// _f_new->Cell(itBound.Down(),8) += _geom->GetInflowVelo()[0];
	// 	// _f_new->Cell(itBound.Down(),7) -= _geom->GetInflowVelo()[0];
	// 	_u->Cell(itBound.Down()) = _geom->GetInflowVelo()[0];
	// }



	/////////////////////////////////////////////////////////////////////////
	////////////// END LATTICE BOLTZMANN IMPLEMENTATION /////////////////////
	/////////////////////////////////////////////////////////////////////////
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
