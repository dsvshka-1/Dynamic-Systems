#pragma once

#include <math.h>
#include <iostream>
#include <fstream>

class DOPRI5
{
    static constexpr double C2 = 0.2; static constexpr double C3 = 0.3; static constexpr double C4 = 0.8; static constexpr double C5 = 8.0 / 9.0; static constexpr double A21 = 0.2; static constexpr double A31 = 3.0 / 40.0; static constexpr double A32 = 9.0 / 40.0;
    static constexpr double A41 = 44.0 / 45.0; static constexpr double A42 = -56.0 / 15.0; static constexpr double A43 = 32.0 / 9.0; static constexpr double A51 = 19372.0 / 6561.0; static constexpr double A52 = -25360.0 / 2187.0;
    static constexpr double A53 = 64448.0 / 6561.0; static constexpr double A54 = -212.0 / 729.0; static constexpr double A61 = 9017.0 / 3168.0; static constexpr double A62 = -355.0 / 33.0; static constexpr double A63 = 46732.0 / 5247.0;
    static constexpr double A64 = 49.0 / 176.0; static constexpr double A65 = -5103.0 / 18656.0; static constexpr double A71 = 35.0 / 384.0; static constexpr double A73 = 500.0 / 1113.0; static constexpr double A74 = 125.0 / 192.0;
    static constexpr double A75 = -2187.0 / 6784.0; static constexpr double A76 = 11.0 / 84.0; static constexpr double E1 = 71.0 / 57600.0; static constexpr double E3 = -71.0 / 16695.0; static constexpr double E4 = 71.0 / 1920.0;
    static constexpr double E5 = -17253.0 / 339200.0; static constexpr double E6 = 22.0 / 525.0; static constexpr double E7 = -1.0 / 40.0;    
    
    double eps = 1e-5, er = 0.0, fac = 1.0, denom = 0.0, h = 0.1, hmax = 0.1, hnew = 0, tph = 0, t = 0, tend = 0;
    int n = 0, reject = 0;
    double* k1 = 0, * k2 = 0, * k3 = 0, * k4 = 0, * k5 = 0, * x = 0, * y = 0, * p = 0; 

    void (*f)(double, const double*, double*, const double*) = 0;    
    void (*procedure)(DOPRI5*) = 0;
    bool (*repeat)(DOPRI5*) = 0;
    std::ofstream file;

    /* Math Functions */
    template<class T> 
    inline T max(T a1, T a2) { return a1 > a2 ? a1 : a2; }

    template<class T>
    inline T min(T a1, T a2) { return a1 < a2 ? a1 : a2; }

    template<class T>
    inline T max(T a1, T a2, T a3) { return  (a1 > a2) ? (a1 > a3 ? a1 : a3) : (a2 > a3 ? a2 : a3); }

    template<class T>
    inline T abs(T a) { return a > 0 ? a : -a; }

    template<class T>
    inline T sqr(T a) { return a * a; }

    /* initialization */


    /* repeat */
    static inline bool Tend(DOPRI5* d) { return d->t < d->tend; }

    /* procedure */
    static inline void Empty(DOPRI5* d) {  }    
    
    static inline void SolOut(DOPRI5* d) {
        d->file << d->t;
        for (int i = 0; i < d->n; i++)
            d->file << ' ' << d->x[i];
        d->file << '\n';
    }

public:    
    DOPRI5() { };
    virtual ~DOPRI5() {};   

    void sol();
    void portrait(int n0, int nn0, double t0, const double* x0, const double* p0, void (*f)(double, const double*, double*, const double*), const char* filename0);
};

