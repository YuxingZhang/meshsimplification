#include <vector>
#include <math.h>
#include <assert.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <cstdio>
#include <cstring>

class Vec3f {
	public:

		//Constructors
		Vec3f();
		Vec3f(float x,float y, float z);
		Vec3f(const Vec3f& v);
		//Deconstructor
		virtual ~Vec3f();
	public:
		//Operators

		//Operator []
		inline float& operator [](int index) {
			assert(index>=0&&index<3);
			return _p[index];
		}

		inline const float& operator [](int index) const {
			assert(index>=0&&index<3);
			return _p[index];
		}

		//Operator =
		Vec3f& operator = (const Vec3f& v);

		//Operators +=,-=, *=, /=
		void operator +=(const Vec3f& v);
		void operator +=(float f);
		void operator -=(const Vec3f& v);
		void operator -=(float f);
		void operator *=(const Vec3f& v);
		void operator *=(float f);
		void operator /=(const Vec3f& v);
		void operator /=(float f);

		//Operators +,-.*,/
		Vec3f operator +(const Vec3f&v) const;
		Vec3f operator +(float f) const;
		Vec3f operator -(const Vec3f&v) const;
		Vec3f operator -(float f) const;
		Vec3f operator *(const Vec3f&v) const;
		Vec3f operator *(float f) const;
		Vec3f operator /(const Vec3f&v) const;
		Vec3f operator /(float f) const;

		Vec3f operator -() const;

	public:
		void Normalize();
		float L2Norm_Sqr();
		bool isEqual(Vec3f v);

	public:
		union {
			struct { float _p[3]; };
			struct { float x,y,z; };
			struct { float r,g,b; };
		};
		enum {_len = 3};   
};
