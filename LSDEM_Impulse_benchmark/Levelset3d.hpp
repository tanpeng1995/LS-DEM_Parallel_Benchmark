#ifndef LEVELSET3D_HPP_
#define LEVELSET3D_HPP_
#include "definitions.hpp"
class Levelset3d{
public:
	Levelset3d();
	Levelset3d(const vector<float> & levelset, const size_t & xdim, const size_t & ydim, const size_t & zdim);
	Levelset3d(const Levelset3d& other);
	Levelset3d & operator=(const Levelset3d & other);
	~Levelset3d();
	bool findPenetration(const Vector3d & point, double & value, Vector3d & gradient) const;
	void clearLevelset();
private:
	vector<float> _levelset;
	size_t _xdim;
	size_t _ydim;
	size_t _zdim;
	double getGridValue(size_t & x, size_t & y, size_t & z) const;
};

Levelset3d::Levelset3d() {
	_xdim = 0; _ydim = 0; _zdim = 0;
}
Levelset3d::Levelset3d(const vector<float> & levelset, const size_t & xdim, const size_t & ydim, const size_t & zdim):
	_levelset(levelset), _xdim(xdim), _ydim(ydim), _zdim(zdim) {
	if (_xdim*_ydim*_zdim != _levelset.size()) {
		cout << "ERROR: Level set size not consistent with dimensions" << endl;
	}
}
Levelset3d::Levelset3d(const Levelset3d& other){
	_xdim = other._xdim;
	_ydim = other._ydim;
	_zdim = other._zdim;
	_levelset = other._levelset;
	_levelset.shrink_to_fit();
}
Levelset3d& Levelset3d::operator=(const Levelset3d & other){
	_xdim = other._xdim;
	_ydim = other._ydim;
	_zdim = other._zdim;
	_levelset = other._levelset;
	_levelset.shrink_to_fit();
	return *this;
}
Levelset3d::~Levelset3d(){
	_xdim = 0; _ydim = 0; _zdim = 0;
	_levelset.clear();
	_levelset.shrink_to_fit();
}
// Vector3d & gradient is "normal" in GRAIN3D_H_
bool Levelset3d::findPenetration(const Vector3d & point, double & value, Vector3d & gradient) const {
	double x = point(0); double y = point(1); double z = point(2);
	if (x+1 > (double)_xdim || y+1 > (double)_ydim || z+1 > (double)_zdim || x < 0 || y < 0 || z < 0 )
		return false;
	size_t xr = (size_t)round(x); size_t yr = (size_t)round(y); size_t zr = (size_t)round(z);
	// check if the point is close to the surface, if not, return false
	if (getGridValue(xr,yr,zr) > 1)
		return false;
	// if the point is close to the surface, find the value by performing trilinear interpolation
	size_t x0 	= (size_t)floor(x);
	size_t y0 	= (size_t)floor(y);
	size_t z0 	= (size_t)floor(z);
	size_t x1 	= (size_t)ceil(x);
	size_t y1 	= (size_t)ceil(y);
	size_t z1 	= (size_t)ceil(z);
	double p000 = getGridValue(x0, y0, z0);
	double p001 = getGridValue(x0,y0,z1);
	double p010 = getGridValue(x0, y1, z0);
	double p011 = getGridValue(x0,y1,z1);
	double p101 = getGridValue(x1,y0,z1);

	double xm 	= getGridValue(x1, y0, z0) - p000;
	double ym 	= p010 - p000;
	double zm	  = p001 - p000;
	double xym	= -xm - p010 + getGridValue(x1,y1,z0);
	double xzm	= -xm - p001 + p101;
	double yzm	= -ym - p001 + p011;
	double xyzm = -xym + p001 - p101 - p011 + getGridValue(x1,y1,z1);
	double dx 	= x - double(x0);
	double dy 	= y - double(y0);
	double dz	  = z - double(z0);
	value = p000 + xm*dx + ym*dy + zm*dz + xym*dx*dy + xzm*dx*dz + yzm*dy*dz + xyzm*dx*dy*dz;
	// if the point lies outside the surface, return false
	if (value >= 0)
		return false;
	// finally, if the point lies in the surface, find the gradient, normalize, and return true
	gradient << xm + xym*dy + xzm*dz + xyzm*dy*dz,
					ym + xym*dx + yzm*dz + xyzm*dx*dz,
					zm + xzm*dx + yzm*dy + xyzm*dx*dy;
	gradient /= gradient.norm();
	return true;
}

void Levelset3d::clearLevelset(){
	_levelset.clear();
	_levelset.shrink_to_fit();
}


double Levelset3d::getGridValue(size_t & x, size_t & y, size_t & z) const {
	return (double)_levelset[z*_ydim*_xdim + y*_xdim + x];
}

#endif /* LEVELSET3D_HPP_ */
