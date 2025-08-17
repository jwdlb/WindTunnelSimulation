//
// Created by jackw on 16/08/2025.
//

#ifndef SIMULATIONGPU_H
#define SIMULATIONGPU_H



class simulationGPU {
public:
    simulation(float densityInp, int numXInp, int numYInp, float hInp);

    void simulate(float dt, float gravity, int numIterations);
    void setScene();
    void setObstacle(float xNorm, float yNorm);

    float sampleField(float x, float y, int field); // field: 0=U,1=V,2=M

private:
    inline int gridIndex(int i, int j) const { return i + j * numX; }

    int numX, numY;
    int numCells;
    float h;
    float density;

    std::vector<float> u, v, newU, newV;
    std::vector<float> m, newM;
    std::vector<float> p, s;

    // Core CPU methods (can migrate to GPU later)
    void integrate(float dt, float gravity);
    void solveIncompressability(int numIterations, float dt);
    void advectVel(float dt);
    void advectSmoke(float dt);
    void extrapolate();
    void clearOldPressures();
};



#endif //SIMULATIONGPU_H
