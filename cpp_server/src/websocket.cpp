#include "websocket.hpp"

WebSocketServer::WebSocketServer() {
    m_server.init_asio();
    
    m_server.set_open_handler([this](websocketpp::connection_hdl hdl) {
        std::cout << "New connection established" << std::endl;
    });
    
    m_server.set_close_handler([this](websocketpp::connection_hdl hdl) {
        std::cout << "Connection closed" << std::endl;
    });
}

WebSocketServer::~WebSocketServer() {
    m_server.stop();
}

void WebSocketServer::start() {
    m_server.listen(8000); // port number
    m_server.start_accept();
    
    m_server.run();
}

void WebSocketServer::stop() {
    m_server.stop();
}

void WebSocketServer::send(connection_hdl hdl, const std::string& message) {
    try {
        m_server.send(hdl, message, websocketpp::frame::opcode::text);
    } catch (const websocketpp::lib::error_code& e) {
        std::cout << "Failed to send message: " << e.message() << std::endl;
    }
}
