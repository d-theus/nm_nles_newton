#include<cmath>
#include<vector>
#include<functional>
#include<string>

using std::vector;
using std::function;
using std::string;

namespace space {
	template <int N>
	class V {
		public:
			V (){x = vector<double>(N);}
			V (double *ar) : V(){
				for (int i = 0; i < N; i++) 
				{
					this->x[i] = ar[i];
				}
			}
			V<N> operator* (double k){
				V<N> res;
				res.x = this -> x;
				int i = 0;
				for(double a : this->x){
					res.x[i] *= k;
					++i;
				}
				return res;
			};
			V<N> operator- (){
				return (*this)*(-1);
			};
			V<N> operator+ (V<N> other){
				V<N> res;
				int i;
				for( double a : this -> x){
					res.x[i] = this -> x[i] + other.x[i];
				}
				return res;
			};
			V<N> operator- (V<N> other){
				return (*this) + (-other);
			};
			string ToString(){
				string res;
				char number[20];
				res += "(";
				int i = 0;
				for(double a : this->x){
					sprintf(number, "%lf", a);
					res += string(number);
					if (i < this->x.size()-1)
						res += ", ";
					++i;
				}
				res += ")";
				return res;
			}
			double& At(size_t i){
				return this->x[i];
			}
			vector<double> x;
		private:

	};

	template <int N>
		class F {
			public:
				function<double(V<N>)> f;
				F ();

			private:
				/* data */
		};


} /* space */
