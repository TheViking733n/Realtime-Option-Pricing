#include <iostream>
#include <cmath>

class OptionModel {
public:    
    OptionModel();
    
    ~OptionModel();
    
    double blackScholes(double S, double K, double r, double sigma, double T, char type);
    
    double binomial(double S, double K, double r, double sigma, double T, int steps, char type);
    
    double monteCarlo(double S, double K, double r, double sigma, double T, int simulations, char type);
    
    void calculateGreeks(double S, double K, double r, double sigma, double T, double &delta, double &gamma, double &theta, double &vega, double &rho, char type);
};
