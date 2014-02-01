
#include "../src/FileReader.hh"
#include "../src/FluidSimulator.hh"
#include "../src/ParticleTracer.hh"
#include "../src/StaggeredGrid.hh"
#include "../src/Debug.hh"

#include "math.h"
#include <iostream>

const real dt = 0.04;
const real gx = 0.0;
const real gy = -1.0;

FileReader setupConf()
{
    FileReader conf;
    //////////////////////////////////////////////////////// Breaking dum /////////////////////////////////////////////////////////////
    conf.registerStringParameter("name");
    conf.setParameter("name", "FreeSurfaceTests");

    //////////////////////////////////////////////////////// Initialization ///////////////////////////////////////////////////////
    conf.registerStringParameter("boundary_condition_N");
    conf.registerStringParameter("boundary_condition_S");
    conf.registerStringParameter("boundary_condition_E");
    conf.registerStringParameter("boundary_condition_W");

    conf.registerRealParameter("boundary_velocity_N");
    conf.registerRealParameter("boundary_velocity_S");
    conf.registerRealParameter("boundary_velocity_E");
    conf.registerRealParameter("boundary_velocity_W");

    conf.registerRealParameter("GX");
    conf.setParameter("GX", gx);
    conf.registerRealParameter("GY");
    conf.setParameter("GY", gy);

    conf.registerRealParameter("Re");
    conf.setParameter("Re", 1.0);

    conf.registerRealParameter("U_INIT");
    conf.registerRealParameter("V_INIT");
    conf.registerRealParameter("P_INIT");

    //////////////////////////////////////////////////////// Geometry Data ////////////////////////////////////////////////////////
    conf.registerRealParameter("xlength");
    conf.setParameter("xlength", 3.0);
    conf.registerRealParameter("ylength");
    conf.setParameter("ylength", 3.0);
    conf.registerIntParameter("imax");
    conf.setParameter("imax", 3);
    conf.registerIntParameter("jmax");
    conf.setParameter("jmax", 3);

    conf.registerIntParameter("ppc");

    conf.registerRealParameter("RectangleParticleX1");
    conf.registerRealParameter("RectangleParticleX2");
    conf.registerRealParameter("RectangleParticleY1");
    conf.registerRealParameter("RectangleParticleY2");

    conf.registerRealParameter("CircleParticleX");
    conf.registerRealParameter("CircleParticleY");
    conf.registerRealParameter("CircleParticleR");

    conf.registerRealParameter("RectangleX1");
    conf.registerRealParameter("RectangleY1");
    conf.registerRealParameter("RectangleX2");
    conf.registerRealParameter("RectangleY2");

    conf.registerRealParameter("CircleX");
    conf.registerRealParameter("CircleY");
    conf.registerRealParameter("CircleR");

    //////////////////////////////////////////////////////// Time Data ////////////////////////////////////////////////////////////
    conf.registerRealParameter("dt");
    conf.setParameter("dt", dt);
    conf.registerIntParameter("timesteps");
    conf.registerRealParameter("safetyfactor");

    //////////////////////////////////////////////////////// Pressure Iteration Data //////////////////////////////////////////////
    conf.registerIntParameter("itermax");
    conf.registerRealParameter("eps");
    conf.registerRealParameter("omg");
    conf.registerRealParameter("gamma");
    conf.registerIntParameter("checkfrequency");
    conf.registerIntParameter("normalizationfrequency");
    conf.setParameter("normalizationfrequency", 1);

    //////////////////////////////////////////////////////// VTK Visualization Data ///////////////////////////////////////////////
    conf.registerIntParameter("outputinterval");
    conf.setParameter("outputinterval", 1);

    return conf;
}

FluidSimulator setupSim()
{
    FileReader conf = setupConf();
    return FluidSimulator(conf);
}

void testInit()
{
    setupSim();
    // No CHECK macro necessary we just check if you does not crash
}

void testOneNeighborEast()
{
    real continuity;
    FluidSimulator sim = setupSim();
    sim.tracer().addRectangle(0, 0, 2, 3, 9);

    srand(0);
    for (int i = 1; i <= sim.grid().imax(); ++i)
    {
        for (int j = 1; j < sim.grid().jmax(); ++j)
        {
            sim.grid().u()(i,j) = rand() % 10 + 1;
            sim.grid().v()(i,j) = rand() % 10 + 1;
        }
    }

    sim.tracer().markCells();
    sim.set_UVP_surface(dt, true);

    continuity = (sim.grid().u()(2, 2) - sim.grid().u()(1, 2)) / sim.grid().dx() +
                 (sim.grid().v()(2, 2) - sim.grid().v()(2, 1)) / sim.grid().dy();
    CHECK(fabs(continuity) < 1e-5);
}

void testTwoNeighborSouthEast()
{
    real continuity;
    FluidSimulator sim = setupSim();
    sim.tracer().addRectangle(0, 3, 2, 1, 9);

    srand(0);
    for (int i = 1; i <= sim.grid().imax(); ++i)
    {
        for (int j = 1; j < sim.grid().jmax(); ++j)
        {
            sim.grid().u()(i,j) = rand() % 10 + 1;
            sim.grid().v()(i,j) = rand() % 10 + 1;
        }
    }
    sim.tracer().markCells();
    sim.set_UVP_surface(dt, true);

    continuity = (sim.grid().u()(2, 2) - sim.grid().u()(1, 2)) / sim.grid().dx() +
                 (sim.grid().v()(2, 2) - sim.grid().v()(2, 1)) / sim.grid().dy();
    CHECK(fabs(continuity) < 1e-5);
}

void testTwoNeighborWestEast()
{
    real continuity;
    FluidSimulator sim = setupSim();
    sim.tracer().addRectangle(1, 0, 2, 3, 9);

    srand(0);
    for (int i = 1; i <= sim.grid().imax(); ++i)
    {
        for (int j = 1; j < sim.grid().jmax(); ++j)
        {
            sim.grid().u()(i,j) = rand() % 10 + 1;
            sim.grid().v()(i,j) = rand() % 10 + 1;
        }
    }
    sim.tracer().markCells();

    real u_old1 = sim.grid().u()(2,2);
    real u_old2 = sim.grid().u()(1,2);

    sim.set_UVP_surface(dt, true);

    real dX = dt * gx;
    real dY = dt * gy;

    CHECK(fabs(sim.grid().u()(2,2) - (u_old1 + dX)) < 1e-5);
    CHECK(fabs(sim.grid().u()(1,2) - (u_old2 + dX)) < 1e-5);

    // TODO: why would the continuity not hold here?
    
    // continuity = (sim.grid().u()(2, 2) - sim.grid().u()(1, 2)) / sim.grid().dx() +
    //              (sim.grid().v()(2, 2) - sim.grid().v()(2, 1)) / sim.grid().dy();
    // CHECK(fabs(continuity) < 1e-5);    
}

int main()
{
    testInit();
    std::cout << "[TEST] Init: OK" << std::endl;

    testOneNeighborEast();
    std::cout << "[TEST] One Neighbor East: OK" << std::endl;

    testTwoNeighborSouthEast();
    std::cout << "[TEST] Two Neighbor South East: OK" << std::endl;

    testTwoNeighborWestEast();
    std::cout << "[TEST] Two Neighbor West East: OK" << std::endl;

    return 0;
}
