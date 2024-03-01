#include "option_models.hpp"

OptionModel::OptionModel() {
}

OptionModel::~OptionModel() {
}

// Black-Scholes model
double OptionModel::blackScholes(double S, double K, double r, double sigma, double T, char type) {
    double d1 = (log(S / K) + (r + 0.5 * pow(sigma, 2)) * T) / (sigma * sqrt(T));
    double d2 = d1 - sigma * sqrt(T);

    if (type == 'c') { // Call option
        return S * cdf(d1) - K * exp(-r * T) * cdf(d2);
    } else if (type == 'p') { // Put option
        return K * exp(-r * T) * cdf(-d2) - S * cdf(-d1);
    } else {
        std::cerr << "Invalid option type" << std::endl;
        return 0.0;
    }
}

// Calculate cumulative distribution function (CDF) of standard normal distribution
double OptionModel::cdf(double x) {
    return 0.5 * (1 + erf(x / sqrt(2)));
}

// Binomial model
double OptionModel::binomial(double S, double K, double r, double sigma, double T, int steps, char type) {
    double dt = T / steps;
    double u = exp(sigma * sqrt(dt));
    double d = 1 / u;
    double p = (exp(r * dt) - d) / (u - d);

    double prices[steps + 1];
    prices[0] = S * pow(d, steps);

    for (int i = 1; i <= steps; ++i) {
        prices[i] = prices[i - 1] * u / d;
    }

    double optionPrices[steps + 1];
    for (int i = 0; i <= steps; ++i) {
        if (type == 'C')
            optionPrices[i] = fmax(0, prices[i] - K);
        else
            optionPrices[i] = fmax(0, K - prices[i]);
    }

    for (int j = steps - 1; j >= 0; --j) {
        for (int i = 0; i <= j; ++i) {
            optionPrices[i] = exp(-r * dt) * (p * optionPrices[i + 1] + (1 - p) * optionPrices[i]);
        }
    }

    return optionPrices[0];
}

// Monte Carlo simulation
double OptionModel::monteCarlo(double S, double K, double r, double sigma, double T, int simulations, char type) {
    srand(time(0));
    double dt = T;
    double sum = 0.0;

    for (int i = 0; i < simulations; ++i) {
        double gaussian = rand() / (RAND_MAX + 1.0);
        double random = sqrt(-2.0 * log(gaussian)) * cos(2 * M_PI * rand() / (RAND_MAX + 1.0));
        double ST = S * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * random);

        double payoff = (type == 'C') ? fmax(ST - K, 0) : fmax(K - ST, 0);
        sum += payoff;
    }

    return (sum / simulations) * exp(-r * T);
}

// Calculate option Greeks (Delta, Gamma, Theta, Vega, Rho)
void OptionModel::calculateGreeks(double S, double K, double r, double sigma, double T, double &delta, double &gamma, double &theta, double &vega, double &rho, char type) {
}
