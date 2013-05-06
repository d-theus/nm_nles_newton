#include <iostream>
#include <functional>
#include <vector>
#include <cmath>

using std::function;
using std::vector;
using std::cout;
using std::endl;

class V2 {
public:
	double x,y;
	double norm()const{return sqrt(x*x + y*y);};
	V2(double _x=0, double _y=0) { x = _x; y = _y;};
	V2 operator+ (V2 p2) {return V2(this->x + p2.x, this->y + p2.y);};
	V2 operator- () {return V2(-x,-y);};
	V2 operator- (V2 p2) {return *this +(-p2);};
};
const V2 hx(10e-10,0);
const V2 hy(0,10e-10);

typedef vector<function<double(V2)> > F2;
typedef vector<vector<function<double(V2)> > > MF2;

F2 grad(function<double(V2)> f){
	F2 res(2);
	res[0] = [f](V2 x){
		return (f(x+hx) - f(x))/hx.norm();
	};
	res[1] = [f](V2 x){
		return (f(x+hy) - f(x))/hy.norm();
	};
	return res;
}

MF2 Jacob(function<double(V2)> f1,function<double(V2)> f2){
	MF2 res(2);
	res[0] = grad(f1);
	res[1] = grad(f2);
	return res;
}

double diff(V2 p1, V2 p2) {return fabs(p1.x-p2.x) > fabs(p1.y-p2.y)?fabs(p1.x-p2.x):fabs(p1.y-p2.y);}
double det(const MF2 &m,const V2 &p){
	return m[0][0](p) * m[1][1](p) - m[1][0](p) * m[0][1](p);
}

V2 find_root(F2 f, V2 p0, double eps){
	V2 p1 = p0;
	V2 dp(0,0);
	MF2 J = Jacob(f[0],f[1]);
	MF2 A1(2), A2(2);

	A1[0] = F2(2);
	A1[1] = F2(2);
	A2[0] = F2(2);
	A2[1] = F2(2);

	A1[0][0] = f[0];
	A1[0][1] = J[0][1];
	A2[0][0] = J[0][0];
	A2[0][1] = f[0];
	A1[1][0] = f[1];
	A1[1][1] = J[1][1];
	A2[1][0] = J[1][0];
	A2[1][1] = f[1];
	do{
		p0 = p1;
		p1.x = p0.x - det(A1,p0) / det (J, p0);
		p1.y = p0.y - det(A2,p0) / det (J, p0);
	}while(diff(p0,p1) > eps);
	return p1;
}

int main(int argc, const char *argv[])
{
	F2 f(2);
	f[0] = [](V2 p) {return p.x*p.x + p.y*p.y -1 ;};
	f[1] = [](V2 p) {return p.y - p.x - 0.5 ;};
	V2 result = find_root(f,V2(1,1), 0.00001);
	cout << result.x << ", " << result.y << endl;
	return 0;
}
