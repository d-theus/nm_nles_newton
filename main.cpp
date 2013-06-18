#include <iostream>
#include <functional>
#include <vector>
#include <cmath>
#include<fstream>
#include<limits>
#include<unistd.h>

using std::function;
using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
using std::numeric_limits;

class MF2;
class F2;
class V2;
F2 grad(function<double(V2)>);

class V2 {
public:
	double x,y;
	double norm_eu()const{return sqrt(x*x + y*y);};
	double norm_max()const{return fabs(x)>fabs(y) ? fabs(x) : fabs(y);};
	V2(double _x=0, double _y=0) { x = _x; y = _y;};
	V2(const V2 &v){x = v.x; y = v.y;};
	V2 operator+ (V2 p2)const {return V2(this->x + p2.x, this->y + p2.y);};
	V2 operator- ()const {return V2(-x,-y);};
	V2 operator- (V2 p2)const {return V2(*this) + (-p2);};
	V2 operator* (double k)const{return V2(k*x, k*y);}
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
	function<double(V2)> Det(){
		return [this](V2 arg){
			return
			(this->f11)(arg)*(this->f22)(arg) - (this->f12)(arg)*(this->f21)(arg);
		};
	}
	double Det(V2 v){
		return this->Det()(v);
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
F2 grad(const function<double(V2)> f){
	F2 res(
	[f](V2 x){
		return (f(x+hx) - f(x))/hx.norm_eu();
	},
	[f](V2 x){
		return (f(x+hy) - f(x))/hy.norm_eu();
	});
	return res;
};
double argmin(function<double(double)> f, double eps, double l, double r){
	static double m1,m2;
	while(fabs(r - l) > eps){
		m1 = l + (r-l)/3;
		m2 = r - (r-l)/3;
		if(f(m1) > f(m2))
			l = m2;
		else 
			r = m1;
	}
	return (r+l)/2;
}

template<class VecFun, class VecScal>
V2 gradient_descent(const VecFun f, VecScal p0, double eps,VecScal app_result, double radius, size_t maxItCount){
	function<double(V2)> phi = [f](V2 x){return f.f1(x)*f.f1(x) + f.f2(x)*f.f2(x);};
	V2 p1 = p0;
	int i = 0;

	ofstream ofs;
	ofs.open("xs.dat", ofstream::out);
	ofs << p1.x << "\t" << p1.y << endl;
	double k;
	do {
		p0 = p1;
		k = argmin(
				[phi,p0](double a){return phi(p0 - phi(p0) * a);},
				eps,
				10e-8,
				1 / (grad(phi)(p0)).norm_eu()
			  );
		++i;
		p1 = p0 - grad(phi)(p0)*k;
		ofs << p1.x << "\t" << p1.y << endl;
		if (i > maxItCount || (app_result - p1).norm_eu() > radius){
			cout << "Algorithm does not converge" << endl;
			ofs.close();
			p1 = V2( std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
			ofs.open("xs.dat", ofstream::out);
			break;
		}
	} while ((p0 - p1).norm_max() > eps);
	ofs.close();
	return p1;
}
template<class VecFun, class VecScal>
V2 find_root(VecFun f, VecScal p0, double eps, VecScal app_result, double radius, size_t maxItCount){
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
		p1.x = p0.x - A1.Det(p0) / J.Det(p0);
		p1.y = p0.y - A2.Det(p0) / J.Det(p0);
		ofs << p1.x << "\t" << p1.y << endl;
		if (i > maxItCount || (app_result - p1).norm_eu() > radius){
			cout << "Algorithm does not converge" << endl;
			ofs.close();
			p1 = V2( std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
			ofs.open("xs.dat", ofstream::out);
			break;
		}
	}while((p0-p1).norm_max() > eps);
	ofs.close();
	return p1;
}

int main(int argc, const char *argv[])
{
	F2 f;
	f.f1 = [](V2 p) {return p.x*p.x + p.y*p.y -1 ;};
	f.f2 = [](V2 p) {return p.y - p.x - 0.5 ;};
	V2 app_result(0.4,0.9);
	ofstream of("p0stats.dat");
	double x=-0.1,y=0.5;
	for (;x<0.9;x+=0.1){
		cout << "x = " << x << endl;
		for (;y<1.3;y+=0.1) {
			V2 p00(x,y);
			V2 p0 = gradient_descent<F2,V2>(f, p00 , 0.01, app_result,1.0f, 1000);
			if(numeric_limits<double>::infinity() != p0.x || numeric_limits<double>::infinity() != p0.y){
				V2 result = find_root<F2,V2>(f,p0, 0.00001, app_result, 0.8,1000);
				cout << result.x << ", " << result.y << endl;
				of << endl << x << "\t" << y << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << endl;
				of.close();
				system("wc -l xs.dat | awk '{print $1}' >> p0stats.dat");
				of.open("p0stats.dat",std::ostream::app);
			}
			else {
				cout << "I know it does not converge here" << endl;
				of << endl << x << "\t" << y << "\t" << 0 << "\t" << 0 << "\t" << 0 << "\t" << endl;
				of << 0 << endl;
			}
		}
		y = 0.3;
	}
	of.close();
	double eps_0 = 10e-02;
	double eps = 10e-02;
	V2 p00(0.5,0.8);
	ofstream ofe("epsstats.dat");
	for(;eps > 10e-09;eps /= 10){
		V2 p0 = gradient_descent<F2,V2>(f, p00 , eps_0, app_result,0.8,10e+05);
			if(numeric_limits<double>::infinity() != p0.x || numeric_limits<double>::infinity() != p0.y){
			V2 result = find_root<F2,V2>(f,p0, eps/100, app_result,0.8,10e+06);
		}
		ofe << eps << endl;
		ofe.close();
		system("wc -l xs.dat | awk '{print $1}' >> epsstats.dat");
		ofe.open("epsstats.dat",std::ostream::app);
	}
	ofe.close();
	return 0;
}
