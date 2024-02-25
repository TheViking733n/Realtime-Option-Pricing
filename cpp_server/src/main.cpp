#include <iostream>
#include "websocket.hpp"

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
