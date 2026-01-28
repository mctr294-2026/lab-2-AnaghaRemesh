#include <functional>
#include "roots.hpp"
#include <cmath>

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    double tol = 1e-6;
    const int max_it=1000000;
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0)
    
        return false;
    double c = 0;
    double fc = 0;
    int iter=0;
    while (iter < max_it)
    {
        iter++;
        
        c = (a + b) / 2;
        fc = f(c);
        if (std::abs(fc) <=tol)
        {
            *root = c;
            return true;
        }


        if (std::abs(fc) < 1e-12)
        {
            *root = c;
            return true;
        }
        else if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else if (fb * fc < 0)
        {
            a = c;
            fa = fc;
        }
    }
    *root = c;
    return false;
}

bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{

    const int maxIter = 1000000;
     int iter = 0;
    const double tol = 1e-6;
    double fa = f(a);
    double fb = f(b);
    double c;
    double fc;
    if (fa * fb > 0){
        return false;
    }
    
    while (iter < maxIter)
    {
         c = a - (fa * (b - a)) / (fb - fa);
        fc = f(c);


        if (std::abs(fc) < tol)
        {
            *root = (a + b) / 2;
            return true;
        }

       
        if (std::abs(fc) < tol)
        {

            *root = c;
            return true;
        }
        if (fa * fc < 0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
        iter++;
    }
    
    return false;
}

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    const int maxIterNewton = 1000000;
    const double tol = 1e-6;
    double x0 = c;
   
    for (int i = 0; i < maxIterNewton; i++)
    {

        if (g(x0) == 0){
            return false;
        }

        double x_new = x0 - f(x0) / g(x0);
        if (x_new < a || x_new > b){
            return false;
        }
        if (std::abs(x_new - x0) < tol)
            {
            *root = x_new;
            return true;
        }
    x0 = x_new;
    }
    return false;
}

bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root){
                const int maxIter = 1000000;
                const double tol = 1e-6;

                 double x0 = c;
                double x1 = c + 1;
                double x2;
                for (int i = 0; i < maxIter; i++)
    {
        if ((f(x0)-f(x1)) == 0.0){
            return false;
        }
   
         x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));
      
        if (std::abs(x2 - x1) < tol)
        {
            *root = x2;
            return true;
        }
        x0 = x1;
        x1 = x2;

        
    }
    return false;
}