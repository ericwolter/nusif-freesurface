# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <ctime>
# include <cstring>

using namespace std;

# include "datadef.hpp"
# include "extras.hpp"
# include "init.hpp"
# include "boundary.hpp"
# include "uvp.hpp"
# include "visual.hpp"
# include "surface.hpp"

int main(int argc, char *Inputfile[]);
void timestamp();

//****************************************************************************80

int main(int argc, char *Inputfile[])

//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for NAST2D.
//
//  Discussion:
//
//    The program is normally invoked with one or more input files specified
//    on the command line, as in
//
//      nast2d my_input_file
//
//  Author:
//
//    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer
//
//  Reference:
//
//    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
//    Numerical Simulation in Fluid Dynamics,
//    SIAM, 1998.
//
//  Parameters:
//
//    Commandline, int ARGC, the number of commandline arguments.
//
//    Commandline, char *Inputfile[], a pointer to the commandline arguments.
//    Normally, Inputfile[0] is the name of the program, Inputfile[1] is the
//    name of the first input file, and so on.
//
{
    REAL beta;
    int cycle;
    REAL del_inj;
    REAL del_streak;
    REAL del_trace;
    REAL del_vec;
    REAL delt;
    REAL delx;
    REAL dely;
    REAL eps;
    REAL **F;
    int **FLAG;
    REAL **G;
    REAL gamma;
    REAL GX;
    REAL GY;
    REAL **HEAT;
    int ibound = 0;
    int ifull = 0;
    int imax;
    int imin = 0;
    char infile[30];
    int init_case;
    int isurf = 0;
    int itermax;
    int itersor = 0;
    int jmax;
    int jmin = 0;
    int N;
    REAL omg;
    char outfile[30];
    REAL **P;
    int p_bound;
    struct particleline *Particlelines;
    REAL pos1x;
    REAL pos1y;
    REAL pos2x;
    REAL pos2y;
    int ppc;
    REAL Pr;
    char problem[30];
    REAL **PSI;
    REAL Re;
    REAL res;
    REAL **RHS;
    char streakfile[30];
    REAL t;
    REAL t_end;
    REAL tau;
    REAL **TEMP;
    REAL TI;
    char tracefile[30];
    REAL **U;
    REAL UI;
    REAL **V;
    char vecfile[30];
    REAL VI;
    int wE;
    int wN;
    int write;
    int wS;
    int wW;
    REAL xlength;
    REAL ylength;
    REAL **ZETA;

    timestamp();

    cout << "\n";
    cout << "NAST2D:\n";
    cout << "\n";
    cout << "  Compiled on " << __DATE__ << " at " << __TIME__ << ".\n";
    cout << "\n";
    cout << "  Griebel, Dornseifer and Neunhoeffer code for 2D\n";
    cout << "  Navier Stokes flow.\n";
    //
    //  READ the parameters of the problem.
    //  Stop if problem type or inputfile are not defined.
    //
    if (READ_PARAMETER(Inputfile[1], problem, &xlength, &ylength, &imax,
                       &jmax, &delx, &dely, &t_end, &delt, &tau, &del_trace, &del_inj,
                       &del_streak, &del_vec, vecfile, tracefile, streakfile, infile, outfile,
                       &N, &pos1x, &pos1y, &pos2x, &pos2y, &itermax, &eps, &omg, &gamma,
                       &p_bound, &Re, &Pr, &beta, &GX, &GY, &UI, &VI, &TI, &wW, &wE, &wN, &wS)
            != 0)
    {
        cout << "\n";
        cout << "NAST2D - Fatal error!\n";
        cout << "  Failure in READ_PARAMETER.\n";
        cout << "  Problem type or input file not specified.\n";
        return (1);
    }

    //
    //  Allocate memory for the arrays
    //
    U    = RMATRIX(0, imax + 1, 0, jmax + 1);
    V    = RMATRIX(0, imax + 1, 0, jmax + 1);
    F    = RMATRIX(0, imax + 1, 0, jmax + 1);
    G    = RMATRIX(0, imax + 1, 0, jmax + 1);
    P    = RMATRIX(0, imax + 1, 0, jmax + 1);
    TEMP = RMATRIX(0, imax + 1, 0, jmax + 1);
    PSI  = RMATRIX(0, imax, 0, jmax);
    ZETA = RMATRIX(1, imax - 1, 1, jmax - 1);
    HEAT = RMATRIX(0, imax, 0, jmax);
    RHS  = RMATRIX(0, imax + 1, 0, jmax + 1);
    FLAG = IMATRIX(0, imax + 1, 0, jmax + 1);
    ppc  = 16;
    //
    //  Read initial values from file "infile"
    //
    init_case = READ_bin(U, V, P, TEMP, FLAG, imax, jmax, infile);

    if (0 < init_case)
    {
        cout << "\n";
        cout << "NAST2D - Fatal error!\n";
        cout << "  Error while reading input file.\n";
        return (1);
    }
    //
    //  Set initial values if "infile" was not specified.
    //
    if (init_case < 0)
    {
        INIT_UVP(problem, U, V, P, TEMP, imax, jmax, UI, VI, TI);

        INIT_FLAG(problem, FLAG, imax, jmax, delx, dely, &ibound);
    }
    //
    //  Initialize particles for streaklines or particle tracing
    //
    if (strcmp(streakfile, "none") || strcmp(tracefile, "none"))
    {
        Particlelines = SET_PARTICLES(N, pos1x, pos1y, pos2x, pos2y);
    }
    //
    //  Initialize particles for free boundary problems
    //
    if (!strcmp(problem, "drop") || !strcmp(problem, "dam"))
    {
        Particlelines = INIT_PARTICLES(&N, imax, jmax, delx, dely, ppc,
                                       problem, U, V);
    }

    SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);

    SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);
    //
    //  Time loop
    //
    for (t = 0.0, cycle = 0; t < t_end; t += delt, cycle++)
    {
        // COMP_delt ( &delt, t, imax, jmax, delx, dely, U, V, Re, Pr, tau, &write,
        //   del_trace, del_inj, del_streak, del_vec );
        //
        //  Determine fluid cells for free boundary problems
        //  and set boundary values at free surface.
        //
        if (!strcmp(problem, "drop") ||
                !strcmp(problem, "dam") ||
                !strcmp(problem, "molding") ||
                !strcmp(problem, "wave"))
        {
            MARK_CELLS(FLAG, imax, jmax, delx, dely, &ifull, &isurf,
                       N, Particlelines);

            cout << "FLAG:" << endl;
            for (int j = jmax; j >= 1; --j)
            {
                for (int i = 1; i <= imax; ++i)
                {
                    cout << (FLAG[i][j] & C_E ? 4 : 2) << " ";
                }
                cout << endl;
            }

            SET_UVP_SURFACE(U, V, P, FLAG, GX, GY, imax, jmax, Re, delx,
                            dely, delt);

            cout << "U:" << endl;
            for (int j = jmax; j >= 1; --j)
            {
                for (int i = 1; i <= imax; ++i)
                {
                    cout << U[i][j] << " ";
                }
                cout << endl;
            }
            cout << "V:" << endl;
            for (int j = jmax; j >= 1; --j)
            {
                for (int i = 1; i <= imax; ++i)
                {
                    cout << V[i][j] << " ";
                }
                cout << endl;
            }
            cout << "P:" << endl;
            for (int j = jmax; j >= 1; --j)
            {
                for (int i = 1; i <= imax; ++i)
                {
                    cout << P[i][j] << " ";
                }
                cout << endl;
            }

        }
        else
        {
            ifull = imax * jmax - ibound;
        }
        //
        //  Compute new temperature
        //
        //COMP_TEMP(U, V, TEMP, FLAG, imax, jmax, delt, delx, dely, gamma, Re, Pr);
        //
        //  Compute tentative velocity field (F,G)
        //
        COMP_FG(U, V, TEMP, F, G, FLAG, imax, jmax, delt, delx, dely, GX,
                GY, gamma, Re, beta);
        cout << "F:" << endl;
        for (int j = jmax; j >= 1; --j)
        {
            for (int i = 1; i <= imax; ++i)
            {
                cout << P[i][j] << " ";
            }
            cout << endl;
        }
        cout << "G:" << endl;
        for (int j = jmax; j >= 1; --j)
        {
            for (int i = 1; i <= imax; ++i)
            {
                cout << P[i][j] << " ";
            }
            cout << endl;
        }
        //
        //  Compute right hand side for pressure equation
        //
        cout << "- COMP_RHS" << endl;
        COMP_RHS(F, G, RHS, FLAG, imax, jmax, delt, delx, dely);
        //
        //  Solve the pressure equation by successive over relaxation
        //
        for (int j = jmax; j >= 1; --j)
        {
            for (int i = 1; i <= imax; ++i)
            {
                cout << RHS[i][j] << " ";
            }
            cout << endl;
        }        
        cout << "- POISSON" << endl;
        if (0 < ifull)
        {
            itersor = POISSON(P, RHS, FLAG, imax, jmax, delx, dely,
                              eps, itermax, omg, &res, ifull, p_bound);
        }
        cout << "P:" << endl;
        for (int j = jmax; j >= 1; --j)
        {
            for (int i = 1; i <= imax; ++i)
            {
                cout << P[i][j] << " ";
            }
            cout << endl;
        }        
        cout << "t= "            << t + delt  << "  "
             << "delt= "         << delt    << "  "
             << "iterations = "  << itersor << "  "
             << "res="           << res     << "  "
             << "F,S,B-cells: "  << ifull   << "  "
             << isurf   << "  "
             << ibound  << "\n";
        //
        //  Compute the new velocity field
        //
        ADAP_UV(U, V, F, G, P, FLAG, imax, jmax, delt, delx, dely);
        //
        //  Set boundary conditions
        //
        SETBCOND(U, V, P, TEMP, FLAG, imax, jmax, wW, wE, wN, wS);
        //
        //  Set special boundary conditions.
        //  Overwrite preset default values.
        //
        SETSPECBCOND(problem, U, V, P, TEMP, imax, jmax, UI, VI);

        if (!strcmp(problem, "drop") ||
                !strcmp(problem, "dam") ||
                !strcmp(problem, "molding") ||
                !strcmp(problem, "wave"))
        {
            SET_UVP_SURFACE(U, V, P, FLAG, GX, GY, imax, jmax, Re, delx, dely, delt);
        }
        //
        //  Write data for visualization
        //
        if ((write & 8) && strcmp(vecfile, "none"))
        {
            COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);

            COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);

            OUTPUTVEC_bin(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
                          imax, jmax, vecfile);
        }

        if ((write & 8) && strcmp(outfile, "none"))
        {
            WRITE_bin(U, V, P, TEMP, FLAG, imax, jmax, outfile);
        }

        // if (strcmp(tracefile, "none"))
        // {
            cout << "trace" << endl;
            PARTICLE_TRACING(tracefile, 1, imax, jmax, delx, dely, delt, U, V, FLAG,
                             N, Particlelines, 0);
        // }

        if (strcmp(streakfile, "none"))
        {
            STREAKLINES(streakfile, write, imax, jmax, delx, dely, delt, t,
                        U, V, FLAG, N, Particlelines);
        }
    }

    COMPPSIZETA(U, V, PSI, ZETA, FLAG, imax, jmax, delx, dely);

    if (strcmp(vecfile, "none"))
    {
        COMP_HEAT(U, V, TEMP, HEAT, FLAG, Re, Pr, imax, jmax, delx, dely);

        OUTPUTVEC_bin(U, V, P, TEMP, PSI, ZETA, HEAT, FLAG, xlength, ylength,
                      imax, jmax, vecfile);
    }

    if (strcmp(outfile, "none"))
    {
        WRITE_bin(U, V, P, TEMP, FLAG, imax, jmax, outfile);
    }
    //
    //  Write the stream function to a file.
    //
    output_one_real(PSI, imin, imax, jmin, jmax, "psi_data.txt");
    //
    //  Free memory
    //
    FREE_RMATRIX(U, 0, imax + 1, 0, jmax + 1);
    FREE_RMATRIX(V, 0, imax + 1, 0, jmax + 1);
    FREE_RMATRIX(F, 0, imax + 1, 0, jmax + 1);
    FREE_RMATRIX(G, 0, imax + 1, 0, jmax + 1);
    FREE_RMATRIX(P, 0, imax + 1, 0, jmax + 1);
    FREE_RMATRIX(TEMP, 0, imax + 1, 0, jmax + 1);
    FREE_RMATRIX(PSI, 0, imax, 0, jmax);
    FREE_RMATRIX(ZETA, 1, imax - 1, 1, jmax - 1);
    FREE_RMATRIX(HEAT, 0, imax, 0, jmax);
    FREE_RMATRIX(RHS, 0, imax + 1, 0, jmax + 1);
    FREE_IMATRIX(FLAG, 0, imax + 1, 0, jmax + 1);

    cout << "\n";
    cout << "NAST2D:\n";
    cout << "  Normal end of execution.\n";

    cout << "\n";
    timestamp();

    return (0);
}
//****************************************************************************80

void timestamp(void)

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    size_t len;
    time_t now;

    now = time(NULL);
    tm = localtime(&now);

    len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

    cout << time_buffer << "\n";

    return;
# undef TIME_SIZE
}
