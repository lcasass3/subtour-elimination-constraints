#include "Solver.hpp"

enum SubtourEliminationTechnique
{
    MTZ,
    GAVISH_GRAVES,
    DFJ
};
class MILPSolver : public Solver
{

    IloEnv _env;
    IloModel _model;
    IloCplex _cplex;
    IloArray<IloNumVarArray> _x;
    IloNumVarArray _u;
    SubtourEliminationTechnique _subtourEliminationTechnique;

    void solveMethod(Solution *S);
    void addMTZConstraints();
    void addGavishGravesConstraints();
    void addDFJConstraints();

public:
    MILPSolver(Instance *I, SubtourEliminationTechnique subtourEliminationTechnique);
    ~MILPSolver();

    double gap() { return _cplex.getMIPRelativeGap(); };

    Solution *recoversolution();
};