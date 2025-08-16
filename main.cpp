#include <SFML/Graphics.hpp>
#include "simulation.h"
#include <chrono>
#include <iostream>

int main() {
    const int windowWidth = 2000;
    const int windowHeight = 1100;
    const int gridSize = 8;

    int numX = windowWidth / gridSize - 2;
    int numY = windowHeight / gridSize - 2;

    simulation sim(1000.0, numX, numY, 0.01);
    sim.setScene();

    sf::RenderWindow window(sf::VideoMode(windowWidth, windowHeight), "SFML Simulation - Optimized");

    // OPTIMIZATION 1: Create texture for fast pixel drawing
    sf::Texture texture;
    texture.create(numX, numY);
    sf::Sprite sprite;
    sprite.setTexture(texture);
    sprite.setScale(gridSize, gridSize);

    // OPTIMIZATION 2: Pre-allocate pixel buffer
    std::vector<sf::Uint8> pixels(numX * numY * 4); // RGBA format

    int frameCount = 0;
    auto lastTime = std::chrono::high_resolution_clock::now();

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // OPTIMIZATION 3: Time simulation separately
        auto simStart = std::chrono::high_resolution_clock::now();
        sim.simulate(1.0 / 120.0, 0.0, 40); // Reduce iterations to match JS
        auto simEnd = std::chrono::high_resolution_clock::now();

        // OPTIMIZATION 4: Fast pixel buffer update (single loop)
        auto renderStart = std::chrono::high_resolution_clock::now();
        for (int i = 1; i < sim.numX - 1; i++) {
            for (int j = 1; j < sim.numY - 1; j++) {
                int simIdx = sim.gridHelper(i, j);
                sf::Uint8 gray = static_cast<sf::Uint8>(std::clamp(sim.m[simIdx], 0.0, 1.0) * 255);

                // Direct pixel access - much faster than rectangles
                int pixelIdx = ((j-1) * numX + (i-1)) * 4;
                pixels[pixelIdx] = gray;     // Red
                pixels[pixelIdx + 1] = gray; // Green
                pixels[pixelIdx + 2] = gray; // Blue
                pixels[pixelIdx + 3] = 255;  // Alpha
            }
        }

        // OPTIMIZATION 5: Single texture update + draw call
        texture.update(pixels.data());

        window.clear(sf::Color::White);
        window.draw(sprite);  // Single draw call!
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