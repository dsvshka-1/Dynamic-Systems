#pragma once
/*
class ODE {
	int _dim;
	int _par_dim;	
public:
	ODE(int dim, int par_dim) : _dim(dim), _par_dim(par_dim) {};
	virtual ~ODE() {};

	virtual inline void f(double, const double*, double*, const double*) const = 0;
	virtual inline void jac(double, const double*, double*, const double*) const = 0;

	inline int dim() { return _dim; }
	inline int par_dim() { return _par_dim; }
};
*/

typedef int (dimension)();
typedef void (RHS)(double t, const double* x, double* y, const double* p);

template<class T>
struct ODE {	
	ODE() = delete;
	static constexpr dimension* dimp() { return &T::dim; }
	static constexpr dimension* par_dimp() { return &T::par_dim; }
	static inline RHS* rhsp() { return &T::rhs; }
	static inline RHS* jacp() { return &T::jac; }
};



#include "Linear.h"

