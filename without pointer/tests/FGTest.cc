#include "../src/FluidSimulator.hh"

#include "../src/SORSolver.hh"
#include "../src/StaggeredGrid.hh"
#include "../src/FileReader.hh"
#include "../src/Types.hh"
#include "../src/Debug.hh"
#include "../src/Array.hh"


#include <math.h>
#define _USE_MATH_DEFINES
#include <cmath>
const real pi = M_PI;

#include <iostream>

#include <algorithm>


void initOwnGridSetup(StaggeredGrid &grid)
{
    PROG("initialize ownGridSetup");
    // Setup 1:
    //    - grid.p   : init with random values
    for (int i = 0; i < grid.p().getSize(0); i++)
    {
        for (int j = 0; j < grid.p().getSize(1); j++)
            grid.p()(i, j) = sin(4 * pi * i * grid.dx()) + sin(4 * pi * j * grid.dy());
    }

    //    - grid.rhs : init with zero
    grid.rhs().fill(0);
}

void initGridSetup3(StaggeredGrid &grid)
{
    PROG("initialize GridSetup2");
    // Setup 2:
    //    - grid.p   : init with random values
    for (int i = 0; i < grid.p().getSize(0); i++)
    {
        for (int j = 0; j < grid.p().getSize(1); j++)
            grid.p()(i, j) = sin(4 * pi * i * grid.dx()) + sin(4 * pi * j * grid.dy());
    }
    //    - grid.rhs : f(x,y) = sin(2 * x * \pi)
    for (int i = 0; i < grid.rhs().getSize(0); i++)
    {
        for (int j = 0; j < grid.rhs().getSize(1); j++)
            grid.rhs()(i, j) = sin(2 * i * grid.dx() * pi);
    }

    // initialize u
    for (int i = 0; i < grid.u().getSize(0); i++)
    {
        for (int j = 0; j < grid.u().getSize(1); j++)
            grid.u()(i, j) = i * grid.dx();
    }

    // initialize v
    for (int i = 0; i < grid.v().getSize(0); i++)
    {
        for (int j = 0; j < grid.v().getSize(1); j++)
            grid.v()(i, j) = j * grid.dy();
    }
}

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        std::cerr << "No config file given" << std::endl;
        return EXIT_FAILURE;
    }

    ////////////////////////////////////////////////////////////// FluidSimulator //////////////////////////////////////////////////////////////

    FileReader confi;

    confi.registerStringParameter("name");

    //////////////////////////////////////////////////////// Initialization ///////////////////////////////////////////////////////
    confi.registerStringParameter("boundary_condition_N");
    confi.registerStringParameter("boundary_condition_S");
    confi.registerStringParameter("boundary_condition_E");
    confi.registerStringParameter("boundary_condition_W");

    confi.registerRealParameter("boundary_velocity_N");
    confi.registerRealParameter("boundary_velocity_S");
    confi.registerRealParameter("boundary_velocity_E");
    confi.registerRealParameter("boundary_velocity_W");

    confi.registerRealParameter("GX");
    confi.registerRealParameter("GY");

    confi.registerRealParameter("Re");

    confi.registerRealParameter("U_INIT");
    confi.registerRealParameter("V_INIT");
    confi.registerRealParameter("P_INIT");

    //////////////////////////////////////////////////////// Geometry Data ////////////////////////////////////////////////////////
    confi.registerRealParameter("xlength");
    confi.registerRealParameter("ylength");
    confi.registerIntParameter("imax");
    confi.registerIntParameter("jmax");

    confi.registerRealParameter("RectangleX1");
    confi.registerRealParameter("RectangleY1");
    confi.registerRealParameter("RectangleX2");
    confi.registerRealParameter("RectangleY2");

    confi.registerRealParameter("CircleX");
    confi.registerRealParameter("CircleY");
    confi.registerRealParameter("CircleR");

    //////////////////////////////////////////////////////// Time Data ////////////////////////////////////////////////////////////
    confi.registerRealParameter("dt");
    confi.registerIntParameter("timesteps");
    confi.registerRealParameter("safetyfactor");

    //////////////////////////////////////////////////////// Pressure Iteration Data //////////////////////////////////////////////
    confi.registerIntParameter("itermax");
    confi.registerRealParameter("eps");
    confi.registerRealParameter("omg");
    confi.registerRealParameter("gamma");
    confi.registerIntParameter("checkfrequency");
    confi.registerIntParameter("normalizationfrequency");

    //////////////////////////////////////////////////////// VTK Visualization Data ///////////////////////////////////////////////
    confi.registerIntParameter("outputinterval");

    CHECK_MSG(confi.readFile(argv[1]), "Could not open file " << argv[1] << " which has to be in the current directory.");

    // create simulator
    FluidSimulator simulator(confi);

    initGridSetup3(simulator.grid());

    simulator.testFG();

    std::cout << "\nU = " << std::endl;
    simulator.grid().u().print();
    std::cout << "\nV = " << std::endl;
    simulator.grid().v().print();
    std::cout << "\nF = " << std::endl;
    simulator.grid().f().print();
    std::cout << "\nG = " << std::endl;
    simulator.grid().g().print();

    // computation of F(1,1) for ownDcavity.par:
    // F(1,1) = 0.5 + 0.1*( 0 - 1 + 0.3375 - 0.1125 - 0.5) = 0.3725
    // analogous for F(1,1), G(0,1) and G(1,1)
    std::cout << "\nF = " << std::endl;
    simulator.grid().f().print();
    std::cout << "\nG = " << std::endl;
    simulator.grid().g().print();

    // check some values of F and G for ownDcavity.par
    CHECK(simulator.grid().f()(1, 0) == 0.3725);  // main for-loop (inner F-values)
    CHECK(simulator.grid().f()(0, 0) == 0);  // boundary values for F
    CHECK(simulator.grid().f()(2, 0) == 1);  // boundary values for F
    CHECK(simulator.grid().g()(0, 1) == 0.3725);  // main for-loop (inner G-values)
    CHECK(simulator.grid().g()(0, 0) == 0);  // boundary values for G
    CHECK(simulator.grid().g()(0, 2) == 1);  // boundary values for G

    std::cout << "FluidSimulator passed test for " << argv[1] << std::endl;

    return 0;
}

