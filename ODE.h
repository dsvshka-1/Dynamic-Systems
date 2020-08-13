#pragma once

class ODE
{
private: 
	int eq_dim, par_dim;

public:
	inline virtual int EQ_DIM() { return eq_dim; }
	inline virtual int PAR_DIM() { return par_dim; }

	ODE(int eq_d, int par_d) : eq_dim(eq_d), par_dim(par_d) {}
	virtual ~ODE() {}

	//void (*F)(const double&, const double*&, double*&, const double*&);
	virtual void F(const double&, const double*&, double*&, const double*&) const {};
};

