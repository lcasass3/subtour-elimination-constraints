// Branch and Cut Solver for the TSP with undirected edges
// CPLEX Generic Callbacks implementation (Incumbent)

#pragma once
#include "Solver.hpp"

struct edge
{
    int id;
    int i;
    int j;
    double w;
};

enum SubtourEliminationTechnique
{
    MTZ,
    GAVISH_GRAVES,
    DFJ
};

class BCSolver : public Solver
{

    IloEnv _env;
    IloModel _model;
    IloCplex _cplex;
    IloBoolVarArray _x;
    vector<edge> _edges;
    SubtourEliminationTechnique _subtourEliminationTechnique;

    void solveMethod(Solution *S);

public:
    BCSolver(Instance *I, SubtourEliminationTechnique subtourEliminationTechnique);

    ~BCSolver() { _env.end(); };

    double gap() { return _cplex.getMIPRelativeGap(); }

    Solution *recoversolution();
};

class MyCallBack : public IloCplex::Callback::Function
{

    Instance *_I;
    IloBoolVarArray _x;
    vector<edge> *_edges;
    SubtourEliminationTechnique _subtourEliminationTechnique;

    void addMTZConstraints(const IloCplex::Callback::Context &context);
    void addGavishGravesConstraints(const IloCplex::Callback::Context &context);
    void addDFJConstraints(const IloCplex::Callback::Context &context);

    void addLazyCuts(const IloCplex::Callback::Context &context);

    void searchHeuristicSolution(const IloCplex::Callback::Context &context);

    void ImproveTour(std::vector<int> &tour, double &length);

    int getEdgeId(int i, int j);

public:
    MyCallBack(Instance *I, IloBoolVarArray &x, vector<edge> *E, SubtourEliminationTechnique technique)
        : _I(I), _x(x), _edges(E), _subtourEliminationTechnique(technique) {}

    // The main callback function
    void invoke(const IloCplex::Callback::Context &context);
};