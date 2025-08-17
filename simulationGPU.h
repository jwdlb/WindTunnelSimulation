#ifndef SIMULATIONGPU_H
#define SIMULATIONGPU_H
#include <cuda_runtime.h>
#include <vector>
using namespace std;

class simulationGPU {
public:
    simulationGPU(float densityInp, int numXInp, int numYInp, float hInp);
    ~simulationGPU();

    void simulate(float dt, float gravity, int numIterations);
    void setScene(int shape);
    void runSolveIncompressibility(int numIterations, float dt);
    void runAdvections(float dt);
    void runExtrapolation();
    void runClearOldPressures();
    void runIntegration(float dt, float gravity);
    void getVelocityGrids(vector<float>& u, vector<float>& v);
    void getSmokeDensityGrid(vector<float>& m);
    void getPressureGrid(vector<float>& p);
    void getSolidFluidGrid(vector<float>& s);
    void changeShape(int shape);

    int numX, numY, numCells, numRows;
    float h, density;

    // Device arrays
    float *d_u, *d_v, *d_newU, *d_newV;
    float *d_m, *d_newM;
    float *d_p, *d_s;
};


#endif //SIMULATIONGPU_H
