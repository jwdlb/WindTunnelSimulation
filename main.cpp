#include <SFML/Graphics.hpp>
#include "simulationGPU.h"
#include <chrono>
#include <iostream>
#include <string>

int main() {
    const int windowWidth = 2000;
    const int windowHeight = 1100;
    const int gridSize = 2;

    int numX = windowWidth / gridSize - 2;
    int numY = windowHeight / gridSize - 2;

    int selectedShape = 0;
    simulationGPU sim(1000.0, numX, numY, 0.01);
    sim.setScene(0);

    std::vector<string> shapes = {"Circle", "Ellipse","Square","Wing"};
    std::vector<float> sStore;
    sim.getSolidFluidGrid(sStore);
    bool smokeUI = true;
    bool outline = true;

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "Wind Tunnel Simulation");

    // OPTIMIZATION 1: Create texture for fast pixel drawing
    sf::Texture texture;
    texture.create(numX, numY);
    sf::Sprite sprite;
    sprite.setTexture(texture);
    sprite.setScale(gridSize, gridSize);
    sf::Font font;
    if (!font.loadFromFile("arial.ttf")) {
        std::cerr << "Font failed to load!" << std::endl;
        return -1;
    } else {
        std::cout << "Font loaded successfully" << std::endl;
    }


    // OPTIMIZATION 2: Pre-allocate pixel buffer
    std::vector<sf::Uint8> pixels(numX * numY * 4); // RGBA format

    int frameCount = 0;
    auto lastTime = std::chrono::high_resolution_clock::now();

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::MouseButtonPressed && event.mouseButton.button == sf::Mouse::Left) {
                sf::Vector2i mousePos = sf::Mouse::getPosition(window);

                // Check which button was clicked
                for (int i = 0; i < 4; i++) {
                    sf::FloatRect buttonRect(10 + i * 230, 10, 220, 60);
                    if (buttonRect.contains(static_cast<float>(mousePos.x), static_cast<float>(mousePos.y))) {
                        selectedShape = i;
                        sim.setScene(i);
                        sim.getSolidFluidGrid(sStore);
                        //std::cout << "Selected Shape: " << i + 1 << std::endl;
                    }
                }
                sf::FloatRect buttonRect(10 + 4 * 230, 10, 220, 60);
                if (buttonRect.contains(static_cast<float>(mousePos.x), static_cast<float>(mousePos.y))) {
                    outline = !outline;
                }
            }
        }

        // OPTIMIZATION 3: Time simulation separately
        auto simStart = std::chrono::high_resolution_clock::now();
        sim.simulate(1.0 / 120.0, 0.0, 40); // Reduce iterations to match JS
        auto simEnd = std::chrono::high_resolution_clock::now();

        std::vector<float> mStore;
        sim.getSmokeDensityGrid(mStore);
        std::vector<float> pStore;
        sim.getPressureGrid(pStore);

        // OPTIMIZATION 4: Fast pixel buffer update (single loop)
        auto renderStart = std::chrono::high_resolution_clock::now();
        for (int i = 1; i < sim.numX - 1; i++) {
            for (int j = 1; j < sim.numY - 1; j++) {
                int simIdx = (i * sim.numY) + j;
                sf::Uint8 grey;
                if (smokeUI) {
                    grey = static_cast<sf::Uint8>(std::clamp(mStore[simIdx], 0.0f, 1.0f) * 255);
                }else {
                    grey = static_cast<sf::Uint8>(std::clamp(pStore[simIdx], 0.0f, 1.0f) * 255);
                }
                int pixelIdx = ((j-1) * numX + (i-1)) * 4;

                if (sStore[simIdx] != 0.0f) {
                    // Inside solid object: keep grey
                    pixels[pixelIdx]     = grey;
                    pixels[pixelIdx + 1] = grey;
                    pixels[pixelIdx + 2] = grey;
                    pixels[pixelIdx + 3] = 255;
                } else{
                    // Empty space: check if next to a solid pixel
                    bool isBorder = false;
                    for (int dx = -1; dx <= 1 && !isBorder; dx++) {
                        for (int dy = -1; dy <= 1 && !isBorder; dy++) {
                            int neighborIdx = ((i + dx) * sim.numY) + (j + dy);
                            if (sStore[neighborIdx] != 0.0f) {
                                isBorder = true;
                            }
                        }
                    }

                    if (isBorder && outline) {
                        // Draw blue outline
                        pixels[pixelIdx]     = 0;   // Red
                        pixels[pixelIdx + 1] = 0;   // Green
                        pixels[pixelIdx + 2] = 50;  // Blue
                        pixels[pixelIdx + 3] = 200; // Alpha
                    } else {
                        // Empty background (fully transparent or white)
                        pixels[pixelIdx]     = 255; // White background
                        pixels[pixelIdx + 1] = 255;
                        pixels[pixelIdx + 2] = 255;
                        pixels[pixelIdx + 3] = 255;
                    }
                }
            }
        }

        // OPTIMIZATION 5: Single texture update + draw call
        texture.update(pixels.data());

        window.clear(sf::Color::White);
        window.draw(sprite);  // Single draw call!

        window.clear(sf::Color::White);

        window.draw(sprite); // simulation

        // // Top bar background
        // sf::RectangleShape topBar(sf::Vector2f(windowWidth, 80));
        // topBar.setFillColor(sf::Color(50, 50, 50));
        // topBar.setPosition(0, 0);
        // window.draw(topBar);

        // Draw shape buttons
        for (int i = 0; i < 4; i++) {
            sf::RectangleShape button(sf::Vector2f(220, 60));
            button.setPosition(10 + i * 230, 10);
            button.setFillColor(i == selectedShape ? sf::Color(150, 150, 250) : sf::Color(100, 100, 100));
            window.draw(button);

            sf::Text text;
            text.setFont(font);
            text.setString(shapes[i]);
            text.setCharacterSize(25);
            text.setFillColor(sf::Color::White);
            text.setPosition(70 + i * 230, 22);
            window.draw(text);
        }

        sf::RectangleShape button(sf::Vector2f(220, 60));
        button.setPosition(10 + 4 * 230, 10);
        button.setFillColor(outline ? sf::Color(0, 150, 0) : sf::Color(100, 100, 100));
        window.draw(button);

        sf::Text text;
        text.setFont(font);
        text.setString("Shape Outline");
        text.setCharacterSize(25);
        text.setFillColor(sf::Color::White);
        text.setPosition(40 + 4 * 230, 22);
        window.draw(text);

        window.display();
        auto renderEnd = std::chrono::high_resolution_clock::now();

        // Performance monitoring
        frameCount++;
        if (frameCount % 60 == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            auto frameTime = std::chrono::duration_cast<std::chrono::milliseconds>(now - lastTime).count();
            auto simTime = std::chrono::duration_cast<std::chrono::microseconds>(simEnd - simStart).count();
            auto renderTime = std::chrono::duration_cast<std::chrono::microseconds>(renderEnd - renderStart).count();

            std::cout << "FPS: " << (60000.0 / frameTime)
                      << " | Sim: " << simTime << "μs"
                      << " | Render: " << renderTime << "μs" << std::endl;
            lastTime = now;
        }
    }

    return 0;
}

// BONUS: If you want even more speed, here's a SIMD-optimized version of the pixel loop:

/*
// OPTIMIZATION 6: SIMD pixel processing (compile with -mavx2)
#include <immintrin.h>

void updatePixelsSIMD(const simulation& sim, std::vector<sf::Uint8>& pixels) {
    const int numX = sim.numX - 2;
    const int numY = sim.numY - 2;

    for (int i = 1; i < sim.numX - 1; i++) {
        for (int j = 1; j < sim.numY - 8; j += 8) { // Process 8 pixels at once
            // Load 8 float values from simulation
            __m256 simValues = _mm256_loadu_ps(&sim.m[sim.gridHelper(i, j)]);

            // Clamp to [0,1] and scale to [0,255]
            __m256 clamped = _mm256_min_ps(_mm256_max_ps(simValues, _mm256_setzero_ps()),
                                          _mm256_set1_ps(1.0f));
            __m256 scaled = _mm256_mul_ps(clamped, _mm256_set1_ps(255.0f));

            // Convert to integers
            __m256i gray = _mm256_cvtps_epi32(scaled);

            // Pack and store (this gets complex for RGBA, but gives huge speedup)
            // ... SIMD pixel packing code ...
        }
    }
}
*/