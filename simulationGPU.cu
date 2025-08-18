#include "cuda_runtime.h"
#include "simulationGPU.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace std;

// Initialisation and destruction of class and it's attributes and methods: -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

simulationGPU::simulationGPU(float densityInp, int numXInp, int numYInp, float hInp) {
    density = densityInp;
    numX = numXInp + 2;
    numY = numYInp + 2;
    numCells = numX * numY;
    h = hInp;
    numRows = numY;
    inletVelocity = 2.0f;
    relativeInletHeight = 0.14f;
    shape = 0;

    // Allocate device memory
    cudaMalloc(&d_u,     numCells * sizeof(float));
    cudaMalloc(&d_v,     numCells * sizeof(float));
    cudaMalloc(&d_newU,  numCells * sizeof(float));
    cudaMalloc(&d_newV,  numCells * sizeof(float));
    cudaMalloc(&d_m,     numCells * sizeof(float));
    cudaMalloc(&d_newM,  numCells * sizeof(float));
    cudaMalloc(&d_p,     numCells * sizeof(float));
    cudaMalloc(&d_s,     numCells * sizeof(float));

    // Initialize with zeros
    std::vector<float> temp0(numCells, 0.0f);
    std::vector<float> temp1(numCells, 1.0f);

    cudaMemcpy(d_u, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_m, temp1.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_p, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_s, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_newU, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_newV, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_newM, temp0.data(), numCells * sizeof(float), cudaMemcpyHostToDevice);
}

simulationGPU::~simulationGPU() {
    cudaFree(d_u); // velocity in the horizontal direction
    cudaFree(d_v); // velocity in the vertical direction
    cudaFree(d_newU);
    cudaFree(d_newV);
    cudaFree(d_m); // density of each cell in the grid
    cudaFree(d_newM);
    cudaFree(d_p); // pressure of each cell i the grid
    cudaFree(d_s); // whether each cell is fluid or solid
}

// Helper functions, all used only in GPU: ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

__device__ inline int gridIndex(int i, int j, int numY) {return ((i * numY) + j);}

__device__ float avgU(int i, int j, float *u, int numY) {
    return (u[gridIndex(i,j - 1,numY)] + u[gridIndex(i,j,numY)] + u[gridIndex(i + 1, j - 1,numY)] + u[gridIndex(i + 1,j,numY)]) * 0.25;
}

__device__ float avgV(int i, int j, float *v, int numY) {
    return (v[gridIndex(i,j - 1,numY)] + v[gridIndex(i,j,numY)] + v[gridIndex(i + 1, j - 1,numY)] + v[gridIndex(i + 1,j,numY)]) * 0.25;
}

__device__ float sampleField(float x, float y, int fieldType, float *u, float *v, float *m, float h, int numX, int numY) {
    float h1 = 1.0 / h;
    float h2 = 0.5 * h;

    // Clamp world coordinates to valid range
    x = max(min(x, numX * h), h);
    y = max(min(y, numY * h), h);

    float dx = 0.0;
    float dy = 0.0;
    const float *f = nullptr;

    switch (fieldType) {
        case 0: f = u; dy = h2; break;
        case 1: f = v; dx = h2; break;
        case 2: f = m; dx = h2; dy = h2; break;
    }

    // Convert world coords -> grid indices
    float x0 = min(floor((x - dx) * h1), numX - 1.0);
    float tx = ((x - dx) - x0 * h) * h1;
    float x1 = min(x0 + 1, numX - 1.0);

    float y0 = min(floor((y - dy) * h1), numY - 1.0);
    float ty = ((y - dy) - y0 * h) * h1;
    float y1 = min(y0 + 1, numY - 1.0);

    float sx = 1.0f - tx;
    float sy = 1.0f - ty;

    // Bilinear interpolation using 1D indexing
    return (sx * sy * f[gridIndex(x0,y0,numY)]) +
           (tx * sy * f[gridIndex(x1,y0,numY)]) +
           (tx * ty * f[gridIndex(x1,y1,numY)]) +
           (sx * ty * f[gridIndex(x0,y1,numY)]);
}

__device__ float nacaThickness(float x, float chord, float t) {
    float X = x / chord;
    return 0.5f * t * chord * (0.2969f*sqrtf(X) - 0.1260f*X - 0.3516f*X*X + 0.2843f*X*X*X - 0.1015f*X*X*X*X);
}

// GPU kernels: -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Kernels for the fluid simulations:

__global__ void integrate(float *v, float *s, int numX, int numY, float gravity, float dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= 1 && i < numX && j >= 1 && j < numY - 1) {
        if (s[gridIndex(i,j,numY)] != 0.0f && s[gridIndex(i,j - 1,numY)] != 0.0f) {
            v[gridIndex(i,j,numY)] += gravity * dt;
        }
    }
}

__global__ void clearOldPressures(int numX, int numY, float *p) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;
    if (i >= numX || j >= numY) return;
    p[gridIndex(i,j,numY)] = 0.0f;
}

__global__ void solveIncompressibility(int grid, float cp, float *s, float *u, float *v, float *p,int numX, int numY) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if ((i + j) % 2 != grid) {return;}

    if (i > 0 && i < numX - 1 && j > 0 && j < numY - 1) {
        if (s[gridIndex(i,j,numY)] == 0.0f) {return;}

        float sTemp = s[gridIndex(i,j,numY)];
        float sx0 = s[gridIndex(i - 1,j,numY)];
        float sx1 = s[gridIndex(i + 1,j,numY)];
        float sy0 = s[gridIndex(i,j - 1,numY)];
        float sy1 = s[gridIndex(i,j + 1,numY)];
        sTemp = sx0 + sx1 + sy0 + sy1;
        if (sTemp == 0.0f) {return;}

        float divergence = u[gridIndex(i + 1,j,numY)] - u[gridIndex(i,j,numY)] + v[gridIndex(i,j + 1,numY)] - v[gridIndex(i,j,numY)];
        float dTemp = -divergence / sTemp;
        dTemp *= 1.9;

        p[gridIndex(i,j,numY)] += cp * dTemp;
        u[gridIndex(i,j,numY)] -= sx0 * dTemp;
        u[gridIndex(i + 1,j,numY)] += sx1 * dTemp;
        v[gridIndex(i,j,numY)] -= sy0 * dTemp;
        v[gridIndex(i,j + 1,numY)] += sy1 * dTemp;
    }
}

__global__ void extrapolateU(int numX, int numY, float *u) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numX) return;

    u[gridIndex(i,0,numY)] = u[gridIndex(i,1,numY)];
    u[gridIndex(i,numY - 1,numY)] = u[gridIndex(i,numY - 2,numY)];
}

__global__ void extrapolateV(int numX, int numY, float *v) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= numY) return;

    v[gridIndex(0,i,numY)] = v[gridIndex(1,i,numY)];
    v[gridIndex(numX - 1,i,numY)] = v[gridIndex(numX - 2,i,numY)];
}

__global__ void advectVel(float dt, float h, int numX, int numY, float *u, float *v, float *s, float *m, float *newU, float *newV) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i <= 0 || i >= numX - 1 || j <= 0 || j >= numY - 1) return; // is the thread in the required grid?

    float h2 = 0.5 * h;

    if (s[gridIndex(i,j,numY)] != 0.0f && s[gridIndex(i - 1,j,numY)] != 0.0f) {
        float x = i * h;
        float y = (j * h) + h2;
        float uTemp = u[gridIndex(i,j,numY)];
        float vTemp = avgV(i,j,v,numY);

        x-= (dt * uTemp);
        y-= (dt * vTemp);
        uTemp = sampleField(x,y,0,u,v,m,h,numX,numY);

        newU[gridIndex(i,j,numY)] = uTemp;
    }
    if (s[gridIndex(i,j,numY)] != 0.0f && s[gridIndex(i,j - 1,numY)] != 0.0f) {
        float x = (i * h) + h2;
        float y = j * h;
        float uTemp = avgU(i,j,u,numY);
        float vTemp = v[gridIndex(i,j,numY)];

        x-= (dt * uTemp);
        y-= (dt * vTemp);
        vTemp = sampleField(x,y,1,u,v,m,h,numX,numY);

        newV[gridIndex(i,j,numY)] = vTemp;
    }
}

__global__ void advectSmoke(float dt, float h, float *u, float *v, float *s, float *m, float *newM, int numX, int numY) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i <= 0 || i >= numX - 1 || j <= 0 || j >= numY - 1) return;

    float h2 = 0.5 * h;

    if (s[gridIndex(i,j,numY)] != 0.0f) {
        float uTemp = (u[gridIndex(i,j,numY)] + u[gridIndex(i + 1,j,numY)]) * 0.5;
        float vTemp = (v[gridIndex(i,j,numY)] + v[gridIndex(i,j + 1,numY)]) * 0.5;
        float x = (i * h) + h2 - (dt * uTemp);
        float y = (j * h) + h2 - (dt * vTemp);

        newM[gridIndex(i,j,numY)] = sampleField(x,y,2,u,v,m,h,numX,numY);
    }
}

    // Kernels for setting up the environment in the application:

__global__ void setUpSceneMemory(int numX, int numY, float *s, float *u, float *v, float *m, float inletVelocity, float relativeInletHeight) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= numX || j >= numY) return;

    int idx = gridIndex(i,j,numY);

    float initialVelocity = inletVelocity;

    float sTemp = 1.0f; // Fluid
    if (i == 0 || j == 0 || j == numY -1) {
        sTemp = 0.0; // Solid
    }

    s[idx] = sTemp;

    if (i == 1) {
        u[idx] = initialVelocity;
    }

    // Inlet region marker
    float inletHeight = relativeInletHeight * numY;
    int minHeight = floor((0.5 * numY) - (0.5 * inletHeight));
    int maxHeight = floor((0.5 * numY) + (0.5 * inletHeight));
    if (j >= minHeight && j < maxHeight && i == 0) {
        m[idx] = 0.0f; // mark as smoky for inlet
    }
    if ((j < minHeight || j >= maxHeight) && i == 0) {
        m[idx] = 1.0f;
    }

}

__global__ void setUpCircleObstacle(int numX, int numY, float *s, float *u, float *v, float *m, float xNorm, float yNorm, float h) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i <= 0 || i >= numX - 1 || j <= 0 || j >= numY - 1) return;

    float worldX = xNorm * (numX - 2) * h;
    float worldY = yNorm * (numY - 2) * h;
    float radius = 0.1f * (numX - 2) * h;
    float r2 = radius * radius;

    float cellX = (i + 0.5f) * h;
    float cellY = (j + 0.5f) * h;
    float dx = cellX - worldX;
    float dy = cellY - worldY;

    if ((dx*dx + dy*dy) < (r2)) {
        int idx = gridIndex(i,j,numY);
        s[idx] = 0.0f;   // mark solid
        m[idx] = 1.0f;   // smoke/marker solid

        u[idx] = u[gridIndex(i + 1, j, numY)] = 0.0f;
        v[idx] = v[gridIndex(i, j + 1, numY)] = 0.0f;
    }
}

__global__ void setUpEllipseObstacle(int numX, int numY, float *s, float *u, float *v, float *m, float xNorm, float yNorm, float h) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i <= 0 || i >= numX - 1 || j <= 0 || j >= numY - 1) return;

    float worldX = xNorm * (numX - 2) * h;
    float worldY = yNorm * (numY - 2) * h;
    float radiusX = 0.16f * (numX - 2) * h;
    float radiusY = 0.12f * (numY - 2) * h;

    float cellX = (i + 0.5f) * h;
    float cellY = (j + 0.5f) * h;
    float dx = (cellX - worldX) / radiusX;
    float dy = (cellY - worldY) / radiusY;


    if (dx*dx + dy*dy < 1.0f) {
        int idx = gridIndex(i,j,numY);
        s[idx] = 0.0f;   // mark solid
        m[idx] = 1.0f;   // smoke/marker solid

        u[idx] = u[gridIndex(i + 1, j, numY)] = 0.0f;
        v[idx] = v[gridIndex(i, j + 1, numY)] = 0.0f;
    }
}

__global__ void setUpSquareObstacle(int numX, int numY, float *s, float *u, float *v, float *m, float xNorm, float yNorm, float h) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i <= 0 || i >= numX - 1 || j <= 0 || j >= numY - 1) return;

    // Center of the square in world coordinates
    float worldX = xNorm * (numX - 2) * h;
    float worldY = yNorm * (numY - 2) * h;

    // Half side length in world coordinates
    float halfSide = 0.5f * 0.16f * (numX - 2) * h;

    float cellX = (i + 0.5f) * h;
    float cellY = (j + 0.5f) * h;
    float dx = fabsf(cellX - worldX);
    float dy = fabsf(cellY - worldY);

    // Square test
    if (dx <= halfSide && dy <= halfSide) {
        int idx = gridIndex(i, j, numY);
        s[idx] = 0.0f;   // mark solid
        m[idx] = 1.0f;   // smoke/marker solid

        u[idx] = u[gridIndex(i + 1, j, numY)] = 0.0f;
        v[idx] = v[gridIndex(i, j + 1, numY)] = 0.0f;
    }
}

__global__ void setUpWingObstacle(int numX, int numY, float *s, float *u, float *v, float *m, float xNorm, float yNorm, float h) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i <= 0 || i >= numX - 1 || j <= 0 || j >= numY - 1) return;

    float worldX = xNorm * (numX - 2) * h; // leading edge
    float worldY = yNorm * (numY - 2) * h; // flat bottom reference

    float chord = 0.3f * (numX - 2) * h;      // chord length
    float thickness = 0.5f * (numY - 2) * h; // thick
    float angleDeg = 6.0f;                   // pitch angle
    float angleRad = angleDeg * 3.14159265f / 180.0f;

    // Flip y-axis so j=0 is bottom
    float cellY = (numY - 1 - j + 0.5f) * h;
    float cellX = (i + 0.5f) * h;
    float xLocal = cellX - worldX;

    if (xLocal >= 0.0f && xLocal <= chord) {
        float yBottom = worldY + xLocal * -tanf(angleRad); // tilt the bottom
        float yTop = yBottom + nacaThickness(xLocal, chord, thickness); // curved top

        if (cellY >= yBottom && cellY <= yTop) {
            int idx = gridIndex(i, j, numY);
            s[idx] = 0.0f;   // mark solid
            m[idx] = 1.0f;   // marker solid

            u[idx] = u[gridIndex(i + 1, j, numY)] = 0.0f;
            v[idx] = v[gridIndex(i, j + 1, numY)] = 0.0f;
        }
    }
}

// Runner methods: ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void simulationGPU::runIntegration(float dt, float gravity) {
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((numX + 15) / 16, (numY + 15) / 16);

    integrate<<<numBlocks, threadsPerBlock>>>(d_v, d_s, numX, numY, gravity, dt);
    cudaDeviceSynchronize();
}

void simulationGPU::runClearOldPressures() {
    dim3 threadsPerBlock(16, 16);
    dim3 numBlocks((numX + 15) / 16, (numY + 15) / 16);

    clearOldPressures<<<numBlocks, threadsPerBlock>>>(numX, numY,d_p);
    cudaDeviceSynchronize();
}

void simulationGPU::runSolveIncompressibility(int numIterations, float dt) {
    float cp = density * h / dt;
    dim3 blockSize(16,16);
    dim3 gridSize((numX + 15)/16, (numY + 15)/16);
    for (int iter = 0; iter < numIterations; iter++) {
        solveIncompressibility<<<gridSize, blockSize>>>(0, cp, d_s, d_u, d_v, d_p, numX, numY);
        cudaDeviceSynchronize();

        solveIncompressibility<<<gridSize, blockSize>>>(1, cp, d_s, d_u, d_v, d_p, numX, numY);
        cudaDeviceSynchronize();
    }
}

void simulationGPU::runExtrapolation() {
    dim3 blockSize(128); // 1D thread block
    dim3 numBlocksU((numX + blockSize.x - 1) / blockSize.x);
    dim3 numBlocksV((numY + blockSize.x - 1) / blockSize.x);

    extrapolateU<<<numBlocksU, blockSize>>>(numX, numY, d_u);
    extrapolateV<<<numBlocksV, blockSize>>>(numX, numY, d_v);
    cudaDeviceSynchronize();
}

void simulationGPU::runAdvections(float dt) {
    dim3 blockSize(16,16);
    dim3 gridSize((numX + 15)/16, (numY + 15)/16);

    // --- Copy current fields into "new" buffers at the beginning ---
    cudaMemcpy(d_newU, d_u, numX * numY * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_newV, d_v, numX * numY * sizeof(float), cudaMemcpyDeviceToDevice);
    cudaMemcpy(d_newM, d_m, numX * numY * sizeof(float), cudaMemcpyDeviceToDevice);

    // --- Advection kernels ---
    advectVel<<<gridSize, blockSize>>>(dt, h, numX, numY, d_u, d_v, d_s, d_m, d_newU, d_newV);
    cudaDeviceSynchronize();

    advectSmoke<<<gridSize, blockSize>>>(dt, h, d_u, d_v, d_s, d_m, d_newM, numX, numY);
    cudaDeviceSynchronize();

    std::swap(d_u, d_newU);
    std::swap(d_v, d_newV);
    std::swap(d_m, d_newM);

}


// Getter methods, for UI grids: --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void simulationGPU::getVelocityGrids(vector<float>& u, vector<float>& v) {
    u.resize(numCells);
    v.resize(numCells);

    cudaMemcpy(u.data(), d_u, numCells * sizeof(float), cudaMemcpyDeviceToHost);
    cudaMemcpy(v.data(), d_v, numCells * sizeof(float), cudaMemcpyDeviceToHost);
}

void simulationGPU::getSmokeDensityGrid(vector<float>& m) {
    m.resize(numCells);
    cudaMemcpy(m.data(), d_m, numCells * sizeof(float), cudaMemcpyDeviceToHost);
}

void simulationGPU::getPressureGrid(vector<float>& p) {
    p.resize(numCells);
    cudaMemcpy(p.data(), d_p, numCells * sizeof(float), cudaMemcpyDeviceToHost);
}

void simulationGPU::getSolidFluidGrid(vector<float>& s) {
    s.resize(numCells);
    cudaMemcpy(s.data(), d_s, numCells * sizeof(float), cudaMemcpyDeviceToHost);
}

// Setting up and changing environment / scene: ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void simulationGPU::setScene() {
    dim3 blockSize(16,16);
    dim3 gridSize((numX + 15)/16, (numY + 15)/16);
    // Initialize the scene
    setUpSceneMemory<<<gridSize, blockSize>>>(numX, numY, d_s, d_u, d_v, d_m, inletVelocity, relativeInletHeight);
    cudaDeviceSynchronize();

    // Add obstacle at normalized coordinates (0.3, 0.5)
    switch (shape) {
        case 0:
            setUpCircleObstacle<<<gridSize, blockSize>>>(numX, numY, d_s, d_u, d_v, d_m, 0.3f, 0.5f, h);
            break;
        case 1:
            setUpEllipseObstacle<<<gridSize, blockSize>>>(numX, numY, d_s, d_u, d_v, d_m, 0.3f, 0.5f, h);
            break;
        case 2:
            setUpSquareObstacle<<<gridSize, blockSize>>>(numX, numY, d_s, d_u, d_v, d_m, 0.3f, 0.5f, h);
            break;
        case 3:
            setUpWingObstacle<<<gridSize, blockSize>>>(numX, numY, d_s, d_u, d_v, d_m, 0.3f, 0.5f, h);
            break;
    }
    cudaDeviceSynchronize();
}

void simulationGPU::updateShape(int shapeInp) {
    shape = shapeInp;
    setScene();
}

void simulationGPU::updateInletVel(float inletVelocityInp) {
    inletVelocity = inletVelocityInp;
    setScene();
}

void simulationGPU::updateInletSize(float inletSizeInp) {
    relativeInletHeight = inletSizeInp;
    setScene();
}

// Simulation coordination, singular run: -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void simulationGPU::simulate(float dt, float gravity, int numIterations) {
    //runIntegration(dt, gravity);  don't need gravity for wind tunnel simulations

    runClearOldPressures();

    runSolveIncompressibility(numIterations, dt);

    runExtrapolation();

    runAdvections(dt);

    cudaDeviceSynchronize(); // wait for GPU
}

// Credits:
// Portions of this code were inspired by:
// Copyright 2022 Matthias Müller - Ten Minute Physics
// MIT License
// Website: https://www.matthiasMueller.info/tenMinutePhysics
// YouTube: https://www.youtube.com/c/TenMinutePhysics
