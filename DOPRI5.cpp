#include "DOPRI5.h"

void DOPRI5::sol()
{
    while (repeat(this)) {
        f(t, x, k1, p);
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i] + h * A21 * k1[i];
        }
        f(t + h * C2, y, k2, p);
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i] + h * (A31 * k1[i] + A32 * k2[i]);
        }
        f(t + h * C3, y, k3, p);
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i] + h * (A41 * k1[i] + A42 * k2[i] + A43 * k3[i]);
        }
        f(t + h * C4, y, k4, p);
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i] + h * (A51 * k1[i] + A52 * k2[i] + A53 * k3[i] + A54 * k4[i]);
        }
        f(t + h * C5, y, k5, p);
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i] + h * (A61 * k1[i] + A62 * k2[i] + A63 * k3[i] + A64 * k4[i] + A65 * k5[i]);
        }
        tph = t + h;
        f(tph, y, k2, p);
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i] + h * (A71 * k1[i] + A73 * k3[i] + A74 * k4[i] + A75 * k5[i] + A76 * k2[i]);
        }
        for (int i = 0; i < n; i++)
        {
            k2[i] = E1 * k1[i] + E3 * k3[i] + E4 * k4[i] + E5 * k5[i] + E6 * k2[i];
        }
        f(tph, y, k3, p);
        for (int i = 0; i < n; i++)
        {
            k4[i] = h * (k2[i] + E7 * k3[i]);
        }
        er = 0;
        for (int i = 0; i < n; i++)
        {
            denom = max(1e-5, abs(y[i]), abs(x[i]));
            er += sqr(k4[i] / denom);
        }
        er = sqrt(er / n);
        fac = max(0.2, min(pow(er / eps, 0.2), 5.0));
        hnew = h / fac;
        if (er < eps)
        {
            for (int i = 0; i < n; i++)
            {
                x[i] = y[i];
            }
            t = tph;
            procedure(this);
            if (hnew > hmax) hnew = hmax;
            if (reject) hnew = min(hnew, h);
            reject = 0;
        }
        else
        {
            if (hnew > h) hnew = h;
            //if (double.IsNaN(er)) hnew = 0.6 * h;
            if (reject) hnew *= 0.9;
            reject = 1;
        }
        h = hnew;
    }
}


void DOPRI5::portrait(int n0, int nn0, double t0, const double* x0, const double* p0, void (*f0)(double, const double*, double*, const double*), const char* filename0) {    

    n = n0;
    k1 = new double[n]; k2 = new double[n]; k3 = new double[n]; k4 = new double[n]; k5 = new double[n]; x = new double[n]; y = new double[n];
    p = new double[nn0];

    t = t0;
    for (int i = 0; i < n; i++) {
        x[i] = x0[i];
    }
    for (int i = 0; i < nn0; i++) {
        p[i] = p0[i];
    }
    file.open(filename0);

    file << 't';
    for (int i = 0; i < n; i++)
        file << " x[" << i << ']';
    file << '\n';

    tend = 100;
    f = f0;
    procedure = SolOut;
    repeat = Tend;

    sol();

    file.close();

    delete[] k1; delete[] k2; delete[] k3; delete[] k4; delete[] k5; delete[] x; delete[] y;
    delete[] p;
}

