#include "math.h"
#include <iostream>

#include "ParticleTracer.hh"


ParticleTracer::ParticleTracer ( StaggeredGrid &grid )
    : grid_(grid)
{
}

void ParticleTracer::markCells()
{

}

void ParticleTracer::addRectangle(int x1, int y1, int x2, int y2)
{
    int minX = std::min(x1, x2);
    int maxX = std::max(x1, x2);
    int minY = std::min(y1, y2);
    int maxY = std::max(y1, y2);

    for (int i = minX; i < maxX; ++i)
    {
        for (int j = minY; j < maxY; ++j)
        {
            this->fillCell(i, j, 9);
        }
    }
}

void ParticleTracer::addCircle(int x, int y, int r)
{
    for (int i = -r; i < r; ++i)
    {
        int h = (int)sqrt(r * r - x * x);

        for (int j = -h; j < h; ++j)
        {
            this->fillCell(x + i, y + j, 9);
        }
    }
}

void ParticleTracer::fillCell(int x, int y, int numParticles)
{
    real cellX = x * grid_.dx();
    real cellY = y * grid_.dy();

    int particlesPerSide = (int)(sqrt(numParticles));

    real deltaX = grid_.dx() / (particlesPerSide + 1);
    real deltaY = grid_.dy() / (particlesPerSide + 1);

    for (int i = 1; i <= particlesPerSide; ++i)
    {
        for (int j = 1; j <= particlesPerSide; ++j)
        {
            Particle p(cellX + i * deltaX, cellY + j * deltaY);
            particles_.push_back(p);
        }
    }
}

void ParticleTracer::print()
{
    for (std::vector<Particle>::iterator p = particles_.begin() ; p != particles_.end(); ++p)
    {
        std::cout << "p: [" << p->x() << ", " << p->y() << "]" << std::endl;
    }
}

void ParticleTracer::advanceParticles(real const dt)
{
    for (std::vector<Particle>::iterator p = particles_.begin() ; p != particles_.end(); ++p)
    {
        real u = interpolateU(p->x(), p->y());
        real v = interpolateV(p->x(), p->y());

        p->setX(p->x() + dt * u);
        p->setY(p->y() + dt * v);
    }
}

real ParticleTracer::interpolateU(real x, real y)
{
    return 0.0;
}

real ParticleTracer::interpolateV(real x, real y)
{
    return 0.0;
}


