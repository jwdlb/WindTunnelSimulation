#include "simulation.h"

#include <iostream>
#include <ctime>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

simulation::simulation(double densityInp, int numXInp, int numYInp, double hInp) {
    density = densityInp;
    numX = numXInp + 2;
    numY = numYInp + 2;
    numCells = numX * numY;
    h = hInp;
    numRows = numY;

    u.resize(numCells, 0.0);
    v.resize(numCells, 0.0);
    newU.resize(numCells, 0.0);
    newV.resize(numCells, 0.0);

    p.resize(numCells, 0.0);
    s.resize(numCells, 0.0);

    m.resize(numCells, 1.0);   // m_layer was initialized with 1.0
    newM.resize(numCells, 0.0);
}

inline int simulation::gridHelper(int i,int j) {return ((i * numRows) + j);}

inline double simulation::gridHelperDouble(double i,double j) {return ((i * numRows) + j);}

void simulation::integrate(double dt, double gravity) {
    for (int i = 1; i < numX; i++) {
        for (int j = 1; j < numY - 1 ; j++) {
            if (s[gridHelper(i,j)] != 0.0 && s[gridHelper(i,j - 1)] != 0.0) {
                v[gridHelper(i,j)] += gravity * dt;
            }
        }
    }
}

void simulation::solveIncompressability(int numIterations, double dt) {
    double cp = density * h / dt;
    
    for (int iter = 0; iter < numIterations; iter++) {
        for (int i = 1; i < numX -1; i++) {
            for (int j = 1; j < numY - 1; j++) {
                if (s[gridHelper(i,j)] == 0.0) {
                    continue;
                }
                double sTemp = s[gridHelper(i,j)];
                double sx0 = s[gridHelper(i - 1,j)];
                double sx1 = s[gridHelper(i + 1,j)];
                double sy0 = s[gridHelper(i,j - 1)];
                double sy1 = s[gridHelper(i,j + 1)];
                sTemp = sx0 + sx1 + sy0 + sy1;
                if (sTemp == 0.0) {
                    continue;
                }
                double div = u[gridHelper(i + 1,j)] - u[gridHelper(i,j)] + v[gridHelper(i,j + 1)] - v[gridHelper(i,j)];
                double pTemp = -div / sTemp;
                pTemp *= 1.7;

                p[gridHelper(i,j)] += cp * pTemp;
                u[gridHelper(i,j)] -= sx0 * pTemp;
                u[gridHelper(i + 1,j)] += sx1 * pTemp;
                v[gridHelper(i,j)] -= sy0 * pTemp;
                v[gridHelper(i,j + 1)] += sy1 * pTemp;

            }
        }
    }
}

void simulation::extrapolate() {
    for (int i = 0; i < numX; i++) {
        u[gridHelper(i,0)] = u[gridHelper(i,1)];
        u[gridHelper(i,numY - 1)] = u[gridHelper(i,numY - 2)];
    }
    for (int i = 0; i < numY; i++) {
        u[gridHelper(0,i)] = u[gridHelper(1,i)];
        u[gridHelper(numX - 1,i)] = u[gridHelper(numX - 2,i)];
    }
}

double simulation::sampleField(double x, double y, fieldType field) {
    double h1 = 1.0 / h;
    double h2 = 0.5 * h;

    // Clamp world coordinates to valid range
    x = max(min(x, numX * h), h);
    y = max(min(y, numY * h), h);

    double dx = 0.0;
    double dy = 0.0;
    const vector<double>* f = nullptr;

    switch (field) {
        case U_Field: f = &u; dy = h2; break;
        case V_Field: f = &v; dx = h2; break;
        case S_Field: f = &m; dx = h2; dy = h2; break;
    }

    // Convert world coords -> grid indices
    double x0 = min(floor((x - dx) * h1), numX - 1.0);
    double tx = ((x - dx) - x0 * h) * h1;
    double x1 = min(x0 + 1, numX - 1.0);

    double y0 = min(floor((y - dy) * h1), numY - 1.0);
    double ty = ((y - dy) - y0 * h) * h1;
    double y1 = min(y0 + 1, numY - 1.0);

    double sx = 1.0 - tx;
    double sy = 1.0 - ty;

    // Bilinear interpolation using 1D indexing
    return (sx * sy * (*f)[gridHelper(x0, y0)]) +
           (tx * sy * (*f)[gridHelper(x1, y0)]) +
           (tx * ty * (*f)[gridHelper(x1, y1)]) +
           (sx * ty * (*f)[gridHelper(x0, y1)]);
}


double simulation::avgU(int i, int j) {
    return (u[gridHelper(i,j - 1)] + u[gridHelper(i,j)] + u[gridHelper(i + 1, j - 1)] + u[gridHelper(i + 1,j)]) * 0.25;
}

double simulation::avgV(int i, int j) {
    return (v[gridHelper(i,j - 1)] + v[gridHelper(i,j)] + v[gridHelper(i + 1, j - 1)] + v[gridHelper(i + 1,j)]) * 0.25;
}

void simulation::advectVel(double dt) {
    newU = u;
    newV = v;
    double h2  = 0.5 * h;
    for (int i = 1; i < numX - 1; i++) {
        for (int j = 1; j < numY - 1; j++) {
            // U component:
            if (s[gridHelper(i,j)] != 0.0 && s[gridHelper(i - 1,j)] != 0.0 && j < numY - 1) {
                double x = i * h;
                double y = (j * h) + h2;
                double uTemp = u[gridHelper(i,j)];
                double vTemp = avgV(i, j);

                x -= (dt * uTemp);
                y -= (dt * vTemp);
                uTemp = sampleField(x, y, U_Field);
                newU[gridHelper(i,j)] = uTemp;
            }
            // V component:
            if (s[gridHelper(i,j)] != 0.0 && s[gridHelper(i,j - 1)] != 0.0 && i < numX - 1) {
                double x = (i * h) + h2;
                double y = (j * h);
                double uTemp = avgU(i, j);
                double vTemp = v[gridHelper(i,j)];

                x -= (dt * uTemp);
                y -= (dt * vTemp);
                vTemp = sampleField(x, y, V_Field);
                newV[gridHelper(i,j)] = vTemp;
            }
        }
    }
    u = newU;
    v = newV;
}

void simulation::advectSmoke(double dt) {
    newM = m;
    double h2 = 0.5 * h;

    for (int i = 1; i < numX - 1; i++) {
        for (int j = 1; j < numY - 1; j++) {
            if (s[gridHelper(i,j)] != 0.0) {
                double uTemp = (u[gridHelper(i,j)] + u[gridHelper(i + 1,j)]) * 0.5;
                double vTemp = (v[gridHelper(i,j)] + v[gridHelper(i,j + 1)]) * 0.5;
                double x = (i * h) + h2 - (dt * uTemp);
                double y = (j * h) + h2 - (dt * vTemp);

                newM[gridHelper(i,j)] = sampleField(x, y, S_Field);
            }
        }
    }
    m = newM;
}

void simulation::clearOldPressures() {
    for (int i = 0; i < p.size(); i++) {
        p[i]= 0.0;
    }
}

void simulation::setObstacle(double xNorm, double yNorm) {
    // Convert normalized coordinates (0..1) to world coordinates
    double worldX = xNorm * (numX - 2) * h;
    double worldY = yNorm * (numY - 2) * h;
    double radius = 0.1 * (numX - 2) * h; // radius in world units
    double r2 = radius * radius;
    // Loop over internal grid cells only
    for (int i = 1; i < numX - 1; i++) {
        for (int j = 1; j < numY - 1; j++) {
            double cellX = (i + 0.5) * h;
            double cellY = (j + 0.5) * h;

            double dx = cellX - worldX;
            double dy = cellY - worldY;

            if ((dx*dx + dy*dy) < r2) {
                int idx = gridHelper(i,j);
                // Mark as solid
                s[idx] = 0.0;
                m[idx] = 1.0;

                // Stop velocity inside the obstacle
                u[idx] = u[gridHelper(i + 1,j)] = 0.0;
                v[idx] = v[gridHelper(i,j + 1)] = 0.0;
            }
        }
    }
}


void simulation::setScene() {
    double initialVelocity = 2.0;
    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            double sTemp = 1.0; // Fluid
            if (i == 0 || j == 0 || j == numY -1) {
                sTemp = 0.0; // Solid
            }
            s[gridHelper(i,j)] = sTemp;

            if (i == 1) {
                u[gridHelper(i,j)] = initialVelocity;
            }
        }
    }

    double inletHeight = 0.2 * numY;
    int minHeight = floor((0.5 * numY) - (0.5 * inletHeight));
    int maxHeight = floor((0.5 * numY) + (0.5 * inletHeight));

    for (int i = minHeight; i < maxHeight; i++) {
        m[i] = 0.0;
    }

    setObstacle(0.3,0.5);
}

void simulation::simulate(double dt, double gravity, int numIterations) {
    integrate(dt, gravity);

    clearOldPressures();
    solveIncompressability(numIterations, dt);

    extrapolate();
    advectVel(dt);
    advectSmoke(dt);
}

