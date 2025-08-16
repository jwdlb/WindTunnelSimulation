#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
using namespace std;

enum fieldType {
    U_Field,
    V_Field,
    S_Field
};

class simulation {
public:
    double density;
    int numX;
    int numY;
    int numCells;
    double h;
    vector<vector<double>> u;
    vector<vector<double>> v;
    vector<vector<double>> newU;
    vector<vector<double>> newV;
    vector<vector<double>> p;
    vector<vector<double>> s;
    vector<vector<double>> m;
    vector<vector<double>> newM;

    simulation(double density,int numX,int numY,double h);

    void integrate(double dt,double gravity);

    void solveIncompressability(int numIterations, double dt);

    void extrapolate();

    double sampleField(double x, double y, fieldType field);

    double avgV(int i, int j);

    double avgU(int i, int j);

    void advectVel(double dt);

    void advectSmoke(double dt);

    void clearOldPressures();

    void setScene();

    void setObstacle(double xNorm, double yNorm);

    void simulate(double dt, double gravity, int numIterations);
};



#endif //SIMULATION_H
