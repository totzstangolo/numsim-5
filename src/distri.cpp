#include "distri.hpp"
#include "geometry.hpp"
#include "iterator.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <algorithm>
using namespace std;

/// Constructs a grid based on a geometry
Distri::Distri(const Geometry *geom){

}

/// Constructs a grid based on a geometry with an offset
// @param geom   Geometry information
// @param offset distance of staggered grid point to cell's anchor point
//               (anchor point = lower left corner)
Distri::Distri(const Geometry *geom, const multi_real_t &offset): Distri(geom) {
	_geom = geom;
	_f0 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f1 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f2 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f3 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f4 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f5 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f6 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f7 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	_f8 = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	real_t ex_tmp[9] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
	real_t ey_tmp[9] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
	real_t t_tmp[9] = {4.0/9.0,  1.0/9.0,  1.0/9.0,  1.0/9.0, 1.0/9.0,
    	  1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
	_ex = new real_t[9];
	_ey = new real_t[9];
	_t = new real_t[9];
	_rho = new real_t[(_geom->Size()[0]+2)*(_geom->Size()[1]+2)];
	for(int iCount=0; iCount < 9; iCount++){
		_ex[iCount] = ex_tmp[iCount];
		_ey[iCount] = ey_tmp[iCount];
		_t[iCount]  = t_tmp[iCount];
	}
}

/// Deletes the grid
Distri::~Distri(){
	if(_f1!=nullptr){
		delete[] _f0;
		delete[] _f1;
		delete[] _f2;
		delete[] _f3;
		delete[] _f4;
		delete[] _f5;
		delete[] _f6;
		delete[] _f7;
		delete[] _f8;
		delete[] _ex;
	    delete[] _ey;
	    delete[] _t;
	    delete[] _rho;
	    // delete[] _geom;
	}
}

///     Initializes the grid with a value
void Distri::Initialize(const real_t &rho){
	index_t gridsize = (_geom->Size()[0]+2)*(_geom->Size()[1]+2);
	for(index_t i = 0; i < gridsize; i++){
		_f0[i] = _t[0];
		_f1[i] = _t[1];
		_f2[i] = _t[2];
		_f3[i] = _t[3];
		_f4[i] = _t[4];
		_f5[i] = _t[5];
		_f6[i] = _t[6];
		_f7[i] = _t[7];
		_f8[i] = _t[8];
		_rho[i] = rho;
	}
}

/// Write access to the grid cell at position [it]
real_t &Distri::Cell(const Iterator &it, index_t vel){
	switch(vel){
		case 0:
			return _f0[it];
			break;
		case 1:
			return _f1[it];
			break;
		case 2:
			return _f2[it];
			break;
		case 3:
			return _f3[it];
			break;
		case 4:
			return _f4[it];
			break;
		case 5:
			return _f5[it];
			break;
		case 6:
			return _f6[it];
			break;
		case 7:
			return _f7[it];
			break;
		case 8:
			return _f8[it];
			break;
	}
}
/// Read access to the grid cell at position [it]
const real_t &Distri::Cell(const Iterator &it, index_t vel) const{
	switch(vel){
		case 0:
			return _f0[it];
			break;
		case 1:
			return _f1[it];
			break;
		case 2:
			return _f2[it];
			break;
		case 3:
			return _f3[it];
			break;
		case 4:
			return _f4[it];
			break;
		case 5:
			return _f5[it];
			break;
		case 6:
			return _f6[it];
			break;
		case 7:
			return _f7[it];
			break;
		case 8:
			return _f8[it];
			break;
	}

}

const real_t &Distri::t(index_t vel) const{
	return _t[vel];
}

const real_t &Distri::ex(index_t vel) const{
	return _ex[vel];
}

const real_t &Distri::ey(index_t vel) const{
	return _ey[vel];
}

const real_t &Distri::rho(const Iterator &it) const{
	return _rho[it];
}

real_t &Distri::rho(const Iterator &it){
	return _rho[it];
}

real_t Distri::sum_vel(const Iterator &it){
	real_t res = 0.0;
	for(int iCount=0;iCount<9;iCount++)
		res += Cell(it,iCount);
	return res;
}

real_t Distri::sum_c_vel_x(const Iterator &it){
	real_t res = 0.0;
	for(int iCount=0;iCount<9;iCount++)
		res += _ex[iCount]*Cell(it,iCount);
	return res;
}

real_t Distri::sum_c_vel_y(const Iterator &it){
	real_t res = 0.0;
	for(int iCount=0;iCount<9;iCount++){
		res += _ey[iCount]*Cell(it,iCount);
	}
	return res;
}
