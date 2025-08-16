#include <SFML/Graphics.hpp>
#include "simulation.h"

int main() {
    const int windowWidth = 2000;
    const int windowHeight = 1100;
    const int gridSize = 15;     // Matches simulation h spacing if needed
    const float updateInterval = 0.05f; // 20 FPS for simulation update

    int numX = windowWidth / gridSize - 2; // subtracting 2 for your ghost cells
    int numY = windowHeight / gridSize - 2;

    simulation sim(1000.0, numX, numY, 0.01); // adjust density/h as needed
    sim.setScene();

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "SFML Simulation Grid");


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        sim.simulate(1.0/120.0, 0.0, 100);

        window.clear(sf::Color::White);

        // Draw grid cells based on sim.m
        for (int i = 1; i < sim.numX - 1; i++) {
            for (int j = 1; j < sim.numY - 1; j++) {
                sf::Uint8 gray = static_cast<sf::Uint8>(sim.m[i][j] * 255);
                sf::RectangleShape cell(sf::Vector2f(gridSize, gridSize));
                cell.setPosition((i - 1) * gridSize, (j - 1) * gridSize);
                cell.setFillColor(sf::Color(gray, gray, gray));
                window.draw(cell);
            }
        }

        window.display();
    }

    return 0;
}
