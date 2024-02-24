#include <iostream>
#include "websocket.hpp"

class OptionModel {
public:    
    OptionModel();
    
    ~OptionModel();
    
    double blackScholes(double S, double K, double r, double sigma, double T, char type);
    
    double binomial(double S, double K, double r, double sigma, double T, int steps, char type);
    
    double monteCarlo(double S, double K, double r, double sigma, double T, int simulations, char type);
    
    void calculateGreeks(double S, double K, double r, double sigma, double T, double &delta, double &gamma, double &theta, double &vega, double &rho, char type);
};

int main() {
    WebSocketServer server;

    server.start();

    // Run the server indefinitely (or until stopped)
    std::cout << "WebSocket server started. Listening for connections..." << std::endl;
    while (true) {
        server.m_server.run();
        // This function will block until the server receives a message or connection event.
    }

    return 0;
}