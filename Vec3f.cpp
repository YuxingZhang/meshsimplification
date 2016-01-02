#include "Vec3f.h"

//  Constructors and Deconstructors
Vec3f::Vec3f(void)
{
	memset(_p,0,sizeof(float)*_len);
}

Vec3f::Vec3f(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

Vec3f::Vec3f(const Vec3f &v)
{
	memcpy(_p,v._p,sizeof(float)*_len);
}

Vec3f::~Vec3f(void)
{

}

//////////////////////////////////////////////////////////////////////////
// Operators

Vec3f& Vec3f::operator =( const Vec3f& v)
{
	memcpy(_p,v._p,sizeof(float)*_len);        
	return (*this);
}

void Vec3f::operator +=(const Vec3f& v)
{
	for(int i=0;i<_len;i++)
		_p[i] += v._p[i];
}
void Vec3f::operator +=(float f)
{
	for(int i=0;i<_len;i++)
		_p[i] += f;
}

void Vec3f::operator -=(const Vec3f& v)
{
	for(int i=0;i<_len;i++)
		_p[i] -= v._p[i];
}
void Vec3f::operator -=(float f)
{
	for(int i=0;i<_len;i++)
		_p[i] -= f;
}

void Vec3f::operator *=(const Vec3f& v)
{
	for(int i=0;i<_len;i++)
		_p[i] *= v._p[i];
}
void Vec3f::operator *=(float f)
{
	for(int i=0;i<_len;i++)
		_p[i] *= f;
}

void Vec3f::operator /=(const Vec3f& v)
{
	for(int i=0;i<_len;i++)
		_p[i] /= v._p[i];
}
void Vec3f::operator /=(float f)
{
	for(int i=0;i<_len;i++)
		_p[i] /= f;
}

Vec3f Vec3f::operator +(const Vec3f&v) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] + v[i];
	return res;
}
Vec3f Vec3f::operator +(float f) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] + f;
	return res;
}

Vec3f Vec3f::operator -(const Vec3f&v) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] - v[i];
	return res;
}
Vec3f Vec3f::operator -(float f) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] - f;
	return res;
}

Vec3f Vec3f::operator *(const Vec3f&v) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] * v[i];
	return res;
}
Vec3f Vec3f::operator *(float f) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] * f;
	return res;
}

Vec3f Vec3f::operator /(const Vec3f&v) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] / v[i];
	return res;
}
Vec3f Vec3f::operator /(float f) const
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = (*this)[i] / f;
	return res;
}

Vec3f Vec3f::operator - () const 
{
	Vec3f res;
	for(int i=0;i<_len;i++)
		res[i] = -(*this)[i];
	return res;
}

//////////////////////////////////////////////////////////////////////////
// Other Methods
void Vec3f::Normalize() {
	float fSqr = L2Norm_Sqr();
	if(fSqr>1e-6)
		(*this) *= 1.0f/sqrt(fSqr);
}

bool Vec3f::isEqual(Vec3f v) { 
	if(this->x != v.x)
		return false;
	if(this->y != v.y)
		return false;
	if(this->z != v.z)
		return false;
	return true;
}

float Vec3f::L2Norm_Sqr() {
	return _p[0]*_p[0] + _p[1]*_p[1] + _p[2]*_p[2];
}
