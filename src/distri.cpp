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
	}
}

///     Initializes the grid with a value
void Distri::Initialize(const real_t &value){
	index_t gridsize = (_geom->Size()[0]+2)*(_geom->Size()[1]+2);
	real_t t0 = value*4.0/9.0;
	real_t t1 = value*1.0/9.0;
	real_t t2 = value*1.0/36.0;
	for(index_t i = 0; i < gridsize; i++){
		_f0[i] = t0;
		_f1[i] = t1;
		_f2[i] = t1;
		_f3[i] = t1;
		_f4[i] = t1;
		_f5[i] = t2;
		_f6[i] = t2;
		_f7[i] = t2;
		_f8[i] = t2;
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
