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
    int numRows;
    vector<double> u;
    vector<double> v;
    vector<double> newU;
    vector<double> newV;
    vector<double> p;
    vector<double> s;
    vector<double> m;
    vector<double> newM;

    simulation(double density,int numX,int numY,double h);

    inline int gridHelper(int i,int j);

    inline double gridHelperDouble(double i,double j);

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
