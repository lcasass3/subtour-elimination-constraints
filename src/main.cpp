#include "Instance.hpp"
#include "MILPSolver.hpp"
#include "BCSolver.hpp"

int main(int argc, const char *argv[])
{
    string filename;
    int timlim = 600;
    SubtourEliminationTechnique technique = SubtourEliminationTechnique::MTZ;
    if (argc > 1)
    {
        filename = argv[1];
        if (argc > 2)
        {
            timlim = atoi(argv[2]);
        }
        if (argc > 3)
        {
            string techinqueString = argv[3];
            if (techinqueString == "MTZ")
            {
                technique = SubtourEliminationTechnique::MTZ;
            }
            else if (techinqueString == "GAVISH_GRAVES")
            {
                technique = SubtourEliminationTechnique::GAVISH_GRAVES;
            }
            else if (techinqueString == "DFJ")
            {
                technique = SubtourEliminationTechnique::DFJ;
            }
            else
            {
                cout << "Invalid subtour elimination technique" << endl;
                return 0;
            }
        }
    }
    else
    {
        cout << "Please provide the filename of the instance to solve" << endl;
        return 0;
    }

    Instance *I = NULL;
    try
    {
        I = new Instance(filename);
        // I->print();
    }
    catch (exception &e)
    {
        cout << "Error: " << e.what() << endl;
        return 1;
    }

    // MILP Solver
    // Solver *S = new MILPSolver(I, technique);

    // BC Solver
    Solver *S = new BCSolver(I, technique);

    S->setparam(Solver::TimLim, timlim);

    S->solve();

    Solution *sol = S->recoversolution();

    sol->print();

    sol->isFeasible();

    try
    {
        string tourfile = "Tours/" + I->instancename() + ".opt.tour";
        Solution *sol2 = new Solution(tourfile, I);
        sol2->print();

        cout << "Feasible solution: " << sol->fitness() << endl;
        cout << "Optimal solution: " << sol2->fitness() << endl;
        cout << "Gap: " << 100 * (sol->fitness() - sol2->fitness()) / sol2->fitness() << "%" << endl;

        delete sol2;
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }

    delete S;
    delete sol;
    delete I;

    return 0;
}