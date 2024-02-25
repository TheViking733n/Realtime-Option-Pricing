#include <iostream>
#include <string>
#include <websocketpp/config/asio_no_tls.hpp>
#include <websocketpp/server.hpp>

class WebSocketServer {
public:    
    typedef websocketpp::server<websocketpp::config::asio> server;
    typedef websocketpp::connection_hdl connection_hdl;
    
    WebSocketServer();
    
    ~WebSocketServer();
    
    void start();
    
    void stop();
    
    void send(connection_hdl hdl, const std::string& message);

private:    
    server m_server;
};
