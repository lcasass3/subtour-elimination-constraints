#include "MILPSolver.hpp"

MILPSolver::MILPSolver(Instance *I) : Solver(I, "MILPSolver"), _env(), _model(_env), _cplex(_model), _x(_env), _u(_env), _technique(_subtourEliminationTechnique)
{
    _x = IloArray<IloNumVarArray>(_env, _I->nnodes());
    for (int i = 0; i < _I->nnodes(); i++)
    {
        _x[i] = IloNumVarArray(_env, _I->nnodes(), 0, 1, IloNumVar::Int);
    }
    _u = IloNumVarArray(_env, _I->nnodes(), 2, _I->nnodes(), IloNumVar::Int);

    // Objective function
    IloExpr obj(_env);
    for (int i = 0; i < _I->nnodes(); i++)
    {
        for (int j = 0; j < _I->nnodes(); j++)
        {
            if (i != j)
            {
                obj += _I->travellingtime(_I->Node(i), _I->Node(j)) * _x[i][j];
            }
        }
    }
    _model.add(IloMinimize(_env, obj));
    obj.end();

    // All nodes must be visited
    for (int i = 0; i < _I->nnodes(); i++)
    {
        IloExpr exprin(_env);
        IloExpr exprout(_env);
        for (int j = 0; j < _I->nnodes(); j++)
        {
            if (i != j)
            {
                exprin += _x[i][j];
                exprout += _x[j][i];
            }
        }
        _model.add(IloRange(_env, 1, exprin, 1));
        _model.add(IloRange(_env, 1, exprout, 1));
        exprin.end();
        exprout.end();
    }

    //* Subtour elimination constraints techniques
    switch (_technique)
    {
    case MTZ:
        addMTZConstraints();
        break;
    case GAVISH_GRAVES:
        addGavishGravesConstraints();
        break;
    case DFJ:
        addDFJConstraints();
        break;
    default:
        break;
    }

    _cplex.exportModel("Model.lp");
}

void MILPSolver::addMTZConstraints()
{
    for (int i = 1; i < _I->nnodes(); i++)
    {
        for (int j = 1; j < _I->nnodes(); j++)
        {
            if (i != j)
            {
                _model.add(_u[i] - _u[j] + 1 <= (_I->nnodes() - 1) * (1 - _x[i][j]));
            }
        }
    }
}

void MILPSolver::addGavishGravesConstraints()
{
    IloEnv env = _env;
    IloNumVarArray f(env, _I->nnodes() * _I->nnodes(), 0, IloInfinity, ILOFLOAT);

    for (int i = 1; i < _I->nnodes(); i++)
    {
        IloExpr inflow(env);
        IloExpr outflow(env);
        for (int j = 0; j < _I->nnodes(); j++)
        {
            if (i != j)
            {
                inflow += f[i * _I->nnodes() + j];
                outflow += f[j * _I->nnodes() + i];
            }
        }
        _model.add(inflow - outflow == 1);
        inflow.end();
        outflow.end();
    }

    for (int i = 0; i < _I->nnodes(); i++)
    {
        for (int j = 0; j < _I->nnodes(); j++)
        {
            if (i != j)
            {
                _model.add(f[i * _I->nnodes() + j] <= (_I->nnodes() - 1) * _x[i][j]);
            }
        }
    }

    f.end();
}

void MILPSolver::addDFJConstraints()
{
    IloEnv env = _env;
    IloNumArray Val(env);

    vector<int> rank(_I->nnodes(), -1);

    // We init the ranks
    for (int i = 0; i < _I->nnodes(); i++)
    {
        rank[i] = -(i + 1);
    }

    // We propagate the ranks
    auto it = find_if(rank.begin(), rank.end(), [](int r)
                      { return r < 0; });
    while (it != rank.end())
    {
        int i = std::distance(rank.begin(), it);
        int ri = fabs(rank[i]);
        rank[i] = ri;
        for (int j = 0; j < _I->nnodes(); j++)
        {
            if (i != j)
            {
                int k = getEdgeId(i, j);
                if (Val[k] > 0.5)
                {
                    int rj = fabs(rank[j]);
                    if (rj > ri)
                    {
                        rank[j] = -ri;
                    }
                }
            }
        }
        it = find_if(rank.begin(), rank.end(), [](int r)
                     { return r < 0; });
    }

    // We get the unique ranks
    unordered_set<int> uniqueRanks;
    for (auto r : rank)
        uniqueRanks.insert(r);

    if (uniqueRanks.size() == 1)
        return;

    // We get the subtours for each rank
    IloRangeArray Cuts(env);
    for (auto r : uniqueRanks)
    {
        // We add a subtour elimination constraint
        IloExpr expr(env);
        int size = 0;
        for (int i = 0; i < _I->nnodes(); i++)
        {
            if (rank[i] == r)
            {
                size++;
                for (int j = i + 1; j < _I->nnodes(); j++)
                {
                    if (rank[j] == r)
                    {
                        int k = getEdgeId(i, j);
                        expr += _x[i][j];
                    }
                }
            }
        }
        Cuts.add(expr <= size - 1);
        expr.end();
    }

    // We add the constraints
    if (Cuts.getSize() > 0)
    {
        cout << "Adding " << Cuts.getSize() << " subtour elimination constraints" << endl;
        _model.add(Cuts);
    }
    Cuts.end();

    Val.end();
}

MILPSolver::~MILPSolver()
{
    _env.end();
}

void MILPSolver::solveMethod(Solution *S)
{
    _cplex.setParam(IloCplex::Param::TimeLimit, _timlim);
    _cplex.setParam(IloCplex::Param::Threads, _threads);
    _cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, _mingap);

    if (S != NULL)
    {
        IloNumVarArray Var(_env);
        IloNumArray Val(_env);
        for (auto i : *S)
        {
            for (auto j : *S)
            {
                if (i != j)
                {
                    Var.add(_x[i][j]);
                    Val.add(1);
                }
            }
        }
        _cplex.addMIPStart(Var, Val);
    }

    _cplex.solve();
}

Solution *MILPSolver::recoversolution()
{
    Solution *S = new Solution(_I);
    try
    {
        int i = 0;
        do
        {
            S->push_back(i);
            for (int j = 0; j < _I->nnodes(); j++)
            {
                if (i != j)
                {
                    if (_cplex.getValue(_x[i][j]) > 0.5)
                    {
                        i = j;
                        break;
                    }
                }
            }
        } while (i != 0);
    }
    catch (IloException &e)
    {
        cout << e << endl;
        delete S;
        S = NULL;
    }

    return S;
}