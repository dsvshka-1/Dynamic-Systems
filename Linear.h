#pragma once

#include "ODE.h"

//
//class Linear2D : public ODE {
//public:
//	Linear2D() :ODE(2, 4) {  }
//	virtual ~Linear2D() {  }
//
//	virtual inline void f(double, const double* x, double* y, const double* p) const {
//		y[0] = p[0] * x[0] + p[1] * x[1];
//		y[1] = p[2] * x[0] + p[3] * x[1];
//	}
//
//	virtual inline void jac(double, const double* x, double* y, const double* p) const {
//		y[0] = p[0]; y[1] = p[1];
//		y[2] = p[2]; y[3] = p[3];
//	}
//};

template<int m>
struct Linear : ODE<Linear<m>>
{
	static constexpr int mm = m * m;
	static constexpr int dim() { return m; }
	static constexpr int par_dim() { return mm; }
	static inline void rhs(double, const double* x, double* y, const double* p){
		int k(0);
		for (int i = 0; i < m; i++) {
			y[i] = 0;			
			for (int j = 0; j < m; j++, k++) {
				y[i] += p[k] * x[j];
			}
		}
	}
	static inline void jac(double, const double* x, double* y, const double* p){		
		for (int i = 0; i < mm; i++) {
			y[i] = p[i];
		}
	}
};
