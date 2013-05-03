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
	V2(double _x=0, double _y=0) { x = _x; y = _y;};
	V2 operator- (V2 p2) {return V2(this->x - p2.x, this->y - p2.y);};
};

double diff(V2 p1, V2 p2) {return fabs(p1.x-p2.x) > fabs(p1.y-p2.y)?fabs(p1.x-p2.x):fabs(p1.y-p2.y);}

V2 find_root( 	vector<function<double(V2)> > f,
		vector<vector<function<double(V2)> > > J,
		V2 p0, 
		double eps
	    ){
	V2 p1;
	V2 dp(0,0);
	do{
		p1.x = f[0](p0) + dp.x*J[0][0](p0) + dp.y*J[0][1](p0);
		p1.y = f[0](p0) + dp.x*J[1][0](p0) + dp.y*J[1][1](p0);
		dp = p1 - p0;
	}while(diff(p0,p1) > eps);
	return p1;
}

int main(int argc, const char *argv[])
{
	vector<function<double(V2)> > f(2);
	f[0] = [](V2 p) {return p.x*p.x + p.y*p.y -1 ;};
	f[1] = [](V2 p) {return p.y - p.x - 0.5 ;};

	vector<function<double(V2)> > grad_f0(2);
	grad_f0[0] = [](V2 p) {return 2*p.x;};
	grad_f0[1] = [](V2 p) {return 2*p.y;};
	
	vector<function<double(V2)> > grad_f1(2);
	grad_f1[0] = [](V2 p) {return -1;};
	grad_f1[1] = [](V2 p) {return 1;};

	vector<vector<function<double(V2)> > > J(2);
	J[0] = grad_f0;
	J[1] = grad_f1;

	V2 result = find_root(f,J, V2(1,1), 0.00001);
	cout << result.x << ", " << result.y << endl;

	return 0;
}
