
#include "../src/ParticleTracer.hh"
#include "../src/StaggeredGrid.hh"
#include "../src/Debug.hh"
#include "../src/VTKWriter.hh"

#include "math.h"
#include <iostream>

void testTotalNumber()
{
    StaggeredGrid grid(2, 2, 1, 1);
    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 2.0, 2.0, 0);

    CHECK(tracer.particles().size() == 4 * 9);
}

void testTotalNumberLarge()
{
    StaggeredGrid grid(1.0, 1.0, 1.0 / 30, 1.0 / 30);
    ParticleTracer tracer(&grid);
    tracer.addRectangle(0, 0, 1.0, 1.0, 0);

    CHECK(tracer.particles().size() == 30 * 30 * 9);
}

void testCellCorrect()
{
    StaggeredGrid grid(2, 2, 1, 1);
    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 2.0, 2.0, 0);

    Particle p = tracer.particles()[0];

    CHECK(p.getCellX(grid.dx()) == 1);
    CHECK(p.getCellY(grid.dy()) == 1);
}

void testFirstPositionCorrect()
{
    StaggeredGrid grid(2, 2, 1, 1);
    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 2.0, 2.0, 0);

    Particle p = tracer.particles()[0];
    CHECK(fabs(p.x() - 0.16667) < 1e-5);
    CHECK(fabs(p.y() - 0.16667) < 1e-5);
}

void testCenterPositionCorrect()
{
    StaggeredGrid grid(3, 3, 1, 1);
    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 3.0, 3.0, 0);

    unsigned int centerOffset = 9 * 4 /* jump over cells */ + 4 /* the index of the center particle in a cell */;
    Particle p = tracer.particles()[centerOffset];
    CHECK(fabs(p.x() - 1.5) < 1e-5);
    CHECK(fabs(p.y() - 1.5) < 1e-5);
}

void testCenterZeroInterpolate()
{
    StaggeredGrid grid(3, 3, 1, 1);
    grid.u().fill(0.0);
    grid.v().fill(0.0);

    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 3.0, 3.0, 0);

    unsigned int centerOffset = 9 * 4 /* jump over cells */ + 4 /* the index of the center particle in a cell */;
    Particle p = tracer.particles()[centerOffset];
    real u = tracer.interpolateU(p.x(), p.y());
    real v = tracer.interpolateV(p.x(), p.y());

    CHECK(fabs(u - 0.0) < 1e-5);
    CHECK(fabs(v - 0.0) < 1e-5);
}

void testCenterSingleInterpolate()
{
    StaggeredGrid grid(3, 3, 1, 1);
    grid.u().fill(0.0);
    grid.v().fill(0.0);

    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 3.0, 3.0, 0);

    grid.u()(2, 2) = 1.0;
    grid.v()(2, 2) = 1.0;

    unsigned int centerOffset = 9 * 4 /* jump over cells */ + 4 /* the index of the center particle in a cell */;
    Particle p = tracer.particles()[centerOffset];
    real u = tracer.interpolateU(p.x(), p.y());
    real v = tracer.interpolateV(p.x(), p.y());

    CHECK(fabs(u - 0.5) < 1e-5);
    CHECK(fabs(v - 0.5) < 1e-5);
}

void testInterpolatedCellsCorrect()
{
    StaggeredGrid grid(3, 3, 1, 1);
    grid.u().fill(0.0);
    grid.v().fill(0.0);

    ParticleTracer tracer(&grid);
    tracer.addRectangle(0.0, 0.0, 3.0, 3.0, 0);

    grid.u()(2, 2) = 1.0;
    grid.u()(2, 1) = 1.0;
    grid.u()(2, 3) = 100.0;

    unsigned int centerOffset = 9 * 4 /* jump over cells */ + 3 /* the index of the bottom center particle in a cell */;
    Particle p = tracer.particles()[centerOffset];
    real u = tracer.interpolateU(p.x(), p.y());

    CHECK(fabs(u - 0.5) < 1e-5);
}

void testMarkCellsNoParticles()
{
    StaggeredGrid grid(3, 3, 1, 1);
    grid.u().fill(0.0);
    grid.v().fill(0.0);

    ParticleTracer tracer(&grid);

    tracer.markCells();
    for (int i = 1; i <= 3; ++i)
    {
        for (int j = 1; j <= 3; ++j)
        {
            CHECK(grid.isEmpty(i, j));
        }
    }
}

void testMarkCellsFullParticles()
{
    StaggeredGrid grid(3, 3, 1, 1);
    grid.u().fill(0.0);
    grid.v().fill(0.0);

    ParticleTracer tracer(&grid);
    tracer.addRectangle(0, 0, 3, 3, 0);

    tracer.markCells();
    for (int i = 1; i <= 3; ++i)
    {
        for (int j = 1; j <= 3; ++j)
        {
            CHECK(grid.isFluid(i, j));
        }
    }
}

void testMarkCellsPartialParticles() {

    StaggeredGrid grid(3, 3, 1, 1);
    grid.u().fill(0.0);
    grid.v().fill(0.0);

    ParticleTracer tracer(&grid);

    // leave one column to the right empty
    tracer.addRectangle(0.0, 0.0, 2.0, 3.0, 0);

    tracer.markCells();
    for (int i = 1; i <= 2; ++i)
    {
        for (int j = 1; j <= 3; ++j)
        {
            CHECK(grid.isFluid(i, j));
        }
    }
    for (int i = 3; i <= 3; ++i)
    {
        for (int j = 1; j <= 3; ++j)
        {
            CHECK(grid.isEmpty(i, j));
        }
    }
}

int main()
{
    testTotalNumber();
    std::cout << "[TEST] Total Number Test: OK" << std::endl;

    testTotalNumberLarge();
    std::cout << "[TEST] Total Number Large Test: OK" << std::endl;

    testCellCorrect();
    std::cout << "[TEST] Correct Cell Test: OK" << std::endl;

    testFirstPositionCorrect();
    std::cout << "[TEST] Correct First Position Test: OK" << std::endl;

    testCenterPositionCorrect();
    std::cout << "[TEST] Correct Center Position Test: OK" << std::endl;

    testCenterZeroInterpolate();
    std::cout << "[TEST] Center Zero Interpolation Test: OK" << std::endl;

    testCenterSingleInterpolate();
    std::cout << "[TEST] Center Single Interpolation: OK" << std::endl;

    testInterpolatedCellsCorrect();
    std::cout << "[TEST] Correct Interpolate Cells: OK" << std::endl;

    testMarkCellsNoParticles();
    std::cout << "[TEST] Marked Cells No Particles: OK" << std::endl;

    testMarkCellsFullParticles();
    std::cout << "[TEST] Marked Cells Full Particles: OK" << std::endl;

    testMarkCellsPartialParticles();
    std::cout << "[TEST] Marked Cells Partial Particles: OK" << std::endl;

    return 0;
}
