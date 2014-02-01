
#include "../src/FileReader.hh"
#include "../src/FluidSimulator.hh"
#include "../src/ParticleTracer.hh"
#include "../src/StaggeredGrid.hh"
#include "../src/Debug.hh"

#include "math.h"
#include <iostream>

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
    conf.registerRealParameter("GY");

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
    conf.setParameter("dt", 0.04);
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

FluidSimulator setupSim() {
    FileReader conf = setupConf();
    return FluidSimulator(conf);
}

void testInit()
{
    setupSim();
    // No CHECK macro necessary we just check if you does not crash
}

void testOneNeighborEast() {
    FluidSimulator sim = setupSim();

    sim.grid().u().fill(0.0);
    sim.grid().v().fill(0.0);
    sim.grid().p().fill(0.0);
    sim.grid().rhs().fill(0.0);

    sim.tracer().addRectangle(1,1,2,2,9);
    sim.tracer().markCells();

    sim.grid().obs().print();
}

int main()
{
    testInit();
    std::cout << "[TEST] Init: OK" << std::endl;

    testOneNeighborEast();

    return 0;
}
