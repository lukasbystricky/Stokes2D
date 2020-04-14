/*
 *  Complex.h
 *  
 */

#ifndef _COMPLEX_H
#define _COMPLEX_H

#include <pmmintrin.h>
#include <iostream>
#include <math.h>
/*
#undef _mm_shuffle_pd

#define _mm_shuffle_pd(__A, __B, __C) (static_cast<__m128d>(__builtin_ia32_shufpd (static_cast<__v2df>(__A), static_cast<__v2df>(__B), (__C))))
 */

class Complex {
	friend std::ostream& operator<<(std::ostream& os,const Complex& c);
	friend std::istream& operator>>(std::istream& is,Complex& c);
	friend double real(const Complex& c);
	friend double imag(const Complex& c);
	friend double abs(const Complex& c);
	friend double abs2(const Complex& c);
	friend double arg(const Complex& c);
	friend Complex log(const Complex& c);
	friend Complex exp(const Complex& c);
	friend Complex sin(const Complex& c);
	friend Complex cos(const Complex& c);
    friend Complex tan(const Complex& c);
    friend Complex sinh(const Complex& c);
	friend Complex cosh(const Complex& c);
    //friend Complex tanh(const Complex& c);
    friend Complex asin(const Complex& c);      //Isn't numerically stable
    friend Complex acos(const Complex& c);      //Isn't numerically stable
    friend Complex atan(const Complex& c);
    friend Complex pow(const Complex& c,const Complex& d);
	friend Complex operator-(double d, const Complex& c2);
	friend Complex operator/(double d, const Complex& c2);
public:
	Complex() { z = _mm_setzero_pd(); }
	Complex(double re) { z = _mm_setr_pd(re,0.0); }
	Complex(double re,double im) { z = _mm_setr_pd(re,im); }
	Complex(const Complex& c) { z = c.z; }
	
	Complex& operator=(const Complex& c) { z = c.z; return *this;}
	
	Complex& operator+=(const Complex& c){ z = _mm_add_pd(z,c.z); return *this; }
	Complex& operator+=(double r) { z = _mm_add_pd(z,_mm_setr_pd(r,0.0)); return *this; }
	
	Complex& operator-=(const Complex& c){ z = _mm_sub_pd(z,c.z); return *this; }
	Complex& operator-=(double r) { z = _mm_sub_pd(z,_mm_setr_pd(r,0.0)); return *this; }
	
	Complex& operator*=(const Complex& c){
		double t[2];
		_mm_storeu_pd(t,c.z);
		__m128d tmp = _mm_mul_pd(z,_mm_set1_pd(t[0]));
		z = _mm_shuffle_pd(z,z,1);
		z = _mm_addsub_pd(tmp,_mm_mul_pd(z,_mm_set1_pd(t[1])));
		return *this;	
	}
    
	Complex& operator*=(double r) { z = _mm_mul_pd(z,_mm_set1_pd(r)); return *this; }
	
	Complex& operator/=(const Complex& c) { 
		double t[2];
		__m128d tz = c.z;
		__m128d denom = _mm_mul_pd(tz,tz);
		_mm_storeu_pd(t,z);
		__m128d tmp = _mm_mul_pd(tz,_mm_set1_pd(t[1]));
		denom = _mm_hadd_pd(denom,denom);
		tz = _mm_shuffle_pd(tz, tz, 1);
		tz = _mm_addsub_pd(tmp, _mm_mul_pd(_mm_set1_pd(t[0]), tz));
		tz = _mm_div_pd(_mm_shuffle_pd(tz, tz, 1),denom);
		z = tz;
		return *this; 
	}
	Complex& operator/=(double r) { z = _mm_div_pd(z,_mm_set1_pd(r)); return *this; }
	
	Complex operator-() const {
		double t[2];
		_mm_storeu_pd(t,z);
		return Complex(-t[0],-t[1]);
	}
	
	void conj() { double t[2]; _mm_storeu_pd(t,z); z = _mm_setr_pd(t[0],-t[1]); }
	void sqrt() { 
		double t[2],m;
		_mm_storeu_pd(t,z);
		__m128d tmp = _mm_mul_pd(z,z);
		tmp = _mm_hadd_pd(tmp,tmp);
		_mm_storeh_pd(&m,tmp);
		double a = 0.5*atan2(t[1],t[0]);
		z = _mm_mul_pd(_mm_setr_pd(cos(a),sin(a)),_mm_set1_pd(::exp(0.25*::log(m))));
	}
private:
	__m128d __attribute__((aligned(16))) z;
};
/*Constant complex numbers.*/
const Complex _i(0.0,1.0);

/*----- Polar constructor -----*/
inline Complex polar(double magn, double arg) {
	Complex c(cos(arg),sin(arg));
	return c*=magn;
}

/*----- Binary arithmetic operators -----*/

/*----- Addition -----*/
inline Complex operator+(const Complex& c1,const Complex& c2) {
	return Complex(c1)+=c2;
}
inline Complex operator+(const Complex& c1,double d) {
	return Complex(c1)+=d;
}
inline Complex operator+(double d, const Complex& c2) {
	return Complex(c2)+=d;
}

/*----- Subtraction -----*/
inline Complex operator-(const Complex& c1,const Complex& c2) {
	return Complex(c1)-=c2;
}
inline Complex operator-(const Complex& c1,double d) {
	return Complex(c1)-=d;
}
inline Complex operator-(double d, const Complex& c2) {
	Complex c;
	c.z = _mm_sub_pd(_mm_setr_pd(d,0.0),c2.z); 
	return c;
}

/*----- Multiplication -----*/
inline Complex operator*(const Complex& c1,const Complex& c2) {
	return Complex(c1)*=c2;
}
inline Complex operator*(const Complex& c1,double d) {
	return Complex(c1)*=d;
}
inline Complex operator*(double d, const Complex& c2) {
	return Complex(c2)*=d;
}

/*----- Division -----*/
inline Complex operator/(const Complex& c1,const Complex& c2) {
	return Complex(c1)/=c2;
}
inline Complex operator/(const Complex& c1,double d) {
	return Complex(c1)/=d;
}
inline Complex operator/(double d, const Complex& c2) {
	Complex c(c2);
	c.conj();
	__m128d denom = _mm_mul_pd(c.z,c.z);
	denom = _mm_hadd_pd(denom,denom);
	c.z = _mm_div_pd(c.z,denom);
	return c*=d;
}

/*----- input-output operators -----*/
inline std::ostream& operator<<(std::ostream& os,const Complex& c) {
	double t[2];
	_mm_storeu_pd(t,c.z);
	os << "(" << t[0] << "," << t[1] << ")";
	return os;
}

inline std::istream& operator>>(std::istream& is,Complex& c) {
	double t[2];
	is >> t[0] >> t[1];
	c.z = _mm_loadu_pd(t);
	return is;
}


/*----- Real and imaginary parts -----*/
inline double real(const Complex& c) {
	double t;
	_mm_storel_pd(&t,c.z);
	return t;
}

inline double imag(const Complex& c) {
	double t;
	_mm_storeh_pd(&t,c.z);
	return t;	
}

/*----- Conjugate -----*/
inline Complex conj(const Complex& c) {
	Complex t(c);
	t.conj();
	return t;	
}

/*----- Absolute value -----*/
inline double abs(const Complex& c) {
	double t;
	__m128d tmp = _mm_mul_pd(c.z,c.z);
	tmp = _mm_hadd_pd(tmp,tmp);
	_mm_storeh_pd(&t,tmp);
	return sqrt(t);
}

/*----- Absolute value without square root-----*/
inline double abs2(const Complex& c) {
	double t;
	__m128d tmp = _mm_mul_pd(c.z,c.z);
	tmp = _mm_hadd_pd(tmp,tmp);
	_mm_storeh_pd(&t,tmp);
	return t;
}

/*----- Argument -----*/
inline double arg(const Complex& c) {
	double t[2];
	_mm_storeu_pd(t,c.z);
	return atan2(t[1],t[0]);
}

/*----- Functions -----*/

/*----- Logarithm -----*/
inline Complex log(const Complex& c) {
	double t[2],m;
	_mm_storeu_pd(t,c.z);
	__m128d tmp = _mm_mul_pd(c.z,c.z);
	tmp = _mm_hadd_pd(tmp,tmp);
	_mm_storeh_pd(&m,tmp);
	return Complex(0.5*log(m),atan2(t[1],t[0]));
}

/*----- Exp -----*/
inline Complex exp(const Complex& c) {
	double t[2],m;
	_mm_storeu_pd(t,c.z);
	m = exp(t[0]);
	Complex w(cos(t[1]),sin(t[1]));
	return w*=m;
}

/*----- Sine -----*/
inline Complex sin(const Complex& c) {
	double t[2];
	_mm_storeu_pd(t,c.z);
    Complex c1 = Complex(sin(t[0]),cos(t[0]));
    Complex c2 = Complex(cosh(t[1]),sinh(t[1]));
    c1.z = _mm_mul_pd(c1.z,c2.z);
    return c1;
}

/*----- Cosine -----*/
inline Complex cos(const Complex& c) {
	double t[2];
	_mm_storeu_pd(t,c.z);
    Complex c1 = Complex(cos(t[0]),-sin(t[0]));
    Complex c2 = Complex(cosh(t[1]),sinh(t[1]));
    c1.z = _mm_mul_pd(c1.z,c2.z);
    return c1;
}

/*----- Tan -----*/
inline Complex tan(const Complex& c) {
	double t[2],den;
    _mm_storeu_pd(t,c.z);
    Complex c1 = Complex(cos(t[0]),sinh(t[1]));
    c1.z = _mm_mul_pd(c1.z,c1.z);
    c1.z = _mm_hadd_pd(c1.z,c1.z);
	_mm_storeh_pd(&den,c1.z);
	return Complex(sin(2*t[0]),sinh(2*t[1]))/(2*den);
}

/*----- Hyperbolic sine -----*/
inline Complex sinh(const Complex& c) {
	double t[2];
	_mm_storeu_pd(t,c.z);
    Complex c1 = Complex(cos(t[1]),sin(t[1]));
    Complex c2 = Complex(sinh(t[0]),cosh(t[0]));
    c1.z = _mm_mul_pd(c1.z,c2.z);
    return c1;
}

/*----- Hyperbolic cosine -----*/
inline Complex cosh(const Complex& c) {
	double t[2];
	_mm_storeu_pd(t,c.z);
    Complex c1 = Complex(cos(t[1]),sin(t[1]));
    Complex c2 = Complex(cosh(t[0]),sinh(t[0]));
    c1.z = _mm_mul_pd(c1.z,c2.z);
    return c1;
}

/*----- Powers -----*/
inline Complex pow(const Complex& c,const Complex& d) {
	return exp(d*log(c));
}
inline Complex pow(const Complex& c,double d) {
	return exp(d*log(c));
}
inline Complex pow(double d,const Complex& c) {
	return exp(c*log(d));
}
inline Complex sqrt(const Complex& c) {
	Complex d = Complex(c);
	d.sqrt();
	return d;
}

inline Complex operator^(const Complex& c,const Complex& d) {
	return exp(d*log(c));
}
inline Complex operator^(const Complex& c,double d) {
	return exp(d*log(c));
}
inline Complex operator^(double d,const Complex& c) {
	return exp(c*log(d));
}

/*----- Arcsine -----*/
inline Complex asin(const Complex& c) {
    return -_i*log(_i*c+sqrt(1-c*c));
}

/*----- Arccosine -----*/
inline Complex acos(const Complex& c) {
    return -_i*log(c+_i*sqrt(1-c*c));
}

/*----- Arctan -----*/
inline Complex atan(const Complex& c) {
    return _i/2*(log((1-_i*c)/(1+_i*c)));
}

#endif
