#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include<fstream>

using std::function;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;

class MF2;
class F2;
class V2;
F2 grad(function<double(V2)>);

class V2 {
public:
	double x,y;
	double norm()const{return sqrt(x*x + y*y);};
	V2(double _x=0, double _y=0) { x = _x; y = _y;};
	V2 operator+ (V2 p2) {return V2(this->x + p2.x, this->y + p2.y);};
	V2 operator- () {return V2(-x,-y);};
	V2 operator- (V2 p2) {return *this +(-p2);};
	V2 operator* (double k){return V2(k*x, k*y);}
};
const V2 hx(10e-10,0);
const V2 hy(0,10e-10);

class MF2 {
public:
	function<double(V2)> f11,f12,f21,f22;
	MF2(){}
	MF2 (function<double(V2)> f1,function<double(V2)> f2,
	    function<double(V2)> f3,function<double(V2)> f4){
		this->f11 = f1;
		this->f12 = f2;
		this->f21 = f3;
		this->f22 = f4;
	}
};

class F2 {
public:
	function<double(V2)> f1,f2;
	F2(){}
	F2 (function<double(V2)> f1,function<double(V2)> f2){
		this->f1 = f1;
		this->f2 = f2;
	}
	MF2 jacobian() const {
		MF2 res(
				grad(f1).f1,
				grad(f1).f2,
				grad(f2).f1,
				grad(f2).f2
		       );
		return res;
	}
	V2 operator()(V2 x){V2 res(f1(x),f2(x));return res;};
};
F2 grad(function<double(V2)> f){
	F2 res(
	[f](V2 x){
		return (f(x+hx) - f(x))/hx.norm();
	},
	[f](V2 x){
		return (f(x+hy) - f(x))/hy.norm();
	});
	return res;
};

double diff(V2 p1, V2 p2) {return fabs(p1.x-p2.x) > fabs(p1.y-p2.y)?fabs(p1.x-p2.x):fabs(p1.y-p2.y);}
double det(const MF2 &m,const V2 &p){
	return m.f11(p) * m.f22(p) - m.f12(p) * m.f21(p);
}

V2 gradient_descent(F2 f, V2 p0, double eps){
	function<double(V2)> phi = [f](V2 x){return f.f1(x)*f.f1(x) + f.f2(x)*f.f2(x);};
	V2 p1 = p0;
	double k = 4*eps;
	int i = 0;

	ofstream ofs;
	ofs.open("xs.dat", ofstream::out);
	ofs << p1.x << "\t" << p1.y << endl;
	do {
		++i;
		p0 = p1;
		p1 = p0 - grad(phi)(p0)*k;
		if (phi(p1) <= phi(p0)) 
			k /= 2;
		ofs << p1.x << "\t" << p1.y << endl;
	} while ((p0 - p1).norm() > eps);
	ofs.close();
	cout << "Gradient descent part iteration count: " << i << endl;
	return p1;
}
V2 find_root(F2 f, V2 p0, double eps){
	V2 p1 = p0;
	V2 dp(0,0);
	MF2 J = f.jacobian();
	MF2 A1, A2;

	A1.f11 = f.f1;
	A1.f12 = J.f12;
	A2.f11 = J.f11;
	A2.f12 = f.f1;
	A1.f21 = f.f2;
	A1.f22 = J.f22;
	A2.f21 = J.f21;
	A2.f22 = f.f1;
	int i = 0;
	ofstream ofs;
	ofs.open("xs.dat", ofstream::app);
	do{
		i++;
		p0 = p1;
		p1.x = p0.x - det(A1,p0) / det (J, p0);
		p1.y = p0.y - det(A2,p0) / det (J, p0);
		ofs << p1.x << "\t" << p1.y << endl;
	}while(diff(p0,p1) > eps);
	ofs.close();
	cout << "Newton's part iteration count: " << i << endl;
	return p1;
}

int main(int argc, const char *argv[])
{
	F2 f;
	f.f1 = [](V2 p) {return p.x*p.x + p.y*p.y -1 ;};
	f.f2 = [](V2 p) {return p.y - p.x - 0.5 ;};
	V2 p0 = gradient_descent(f, V2(1,1), 0.01);
	V2 result = find_root(f,p0, 0.00001);
	cout << result.x << ", " << result.y << endl;
	return 0;
}
