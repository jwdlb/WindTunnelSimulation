#include <SFML/Graphics.hpp>
#include "simulation.h"

int main() {
    const int windowWidth = 2000;
    const int windowHeight = 1100;
    const int gridSize = 15;

    int numX = windowWidth / gridSize - 2; // account for ghost cells
    int numY = windowHeight / gridSize - 2;

    simulation sim(1000.0, numX, numY, 0.01);
    sim.setScene();

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "SFML Simulation Grid");

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Advance simulation
        sim.simulate(1.0 / 120.0, 0.0, 100); // dt hard-coded

        window.clear(sf::Color::White);

        // Draw grid cells based on sim.m
        for (int i = 1; i < sim.numX - 1; i++) {
            for (int j = 1; j < sim.numY - 1; j++) {
                int idx = sim.gridHelper(i, j);
                sf::Uint8 gray = static_cast<sf::Uint8>(std::clamp(sim.m[idx], 0.0, 1.0) * 255);

                sf::RectangleShape cell(sf::Vector2f(gridSize, gridSize));
                cell.setPosition(i * gridSize, j * gridSize);
                cell.setFillColor(sf::Color(gray, gray, gray));
                window.draw(cell);
            }
        }

        window.display();
    }

    return 0;
}

