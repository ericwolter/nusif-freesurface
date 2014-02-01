#include "math.h"
#include <iostream>

#include "ParticleTracer.hh"

ParticleTracer::ParticleTracer()
{
}

ParticleTracer::ParticleTracer(StaggeredGrid *grid)
    : grid_(grid)
{
}

void ParticleTracer::markCells()
{
    PROG("marking cells");

    for (int i = 1; i <= grid_->imax(); ++i)
    {
        for (int j = 1; j <= grid_->jmax(); ++j)
        {
            if (grid_->isObstacle(i,j)) continue;

            // grid_->setCellToEmpty(i, j);
        }
    }

    for (std::vector<Particle>::iterator p = particles_.begin() ; p != particles_.end(); ++p)
    {
        int i = p->getCellX(grid_->dx());
        int j = p->getCellY(grid_->dy());

        grid_->setCellToFluid(i, j);
    }

    grid_->refreshEmpty();
}

void ParticleTracer::addRectangle(int x1, int y1, int x2, int y2, int type)
{
    PROG("adding particle rectangle: " << "(" << x1 << "|" << y1 << ")" << ", " << "(" << x2 << "|" << y2 << ")");
    int minX = std::min(x1, x2);
    int maxX = std::max(x1, x2);
    int minY = std::min(y1, y2);
    int maxY = std::max(y1, y2);

    for (int i = minX; i < maxX; ++i)
    {
        for (int j = minY; j < maxY; ++j)
        {
            // std::cout << "TRACER fillCell: " << i + 1 << ", " << j + 1 << std::endl;
            this->fillCell(i + 1, j + 1, grid_->ppc(), type);
        }
    }
}

void ParticleTracer::addCircle(int x, int y, int r, int type)
{
    PROG("adding particle circle: " << "(" << x << "|" << y << "|" << r << ")");

    for (int i = -r; i < r; ++i)
    {
        int h = (int)sqrt(r * r - x * x);

        for (int j = -h; j < h; ++j)
        {
            this->fillCell(x + i, y + j, grid_->ppc(), type);
        }
    }
}

void ParticleTracer::fillCell(int x, int y, int numParticles, int type)
{
    real cellX = (x - 1) * grid_->dx();
    real cellY = (y - 1) * grid_->dy();
    // std::cout << "TRACER cellX,Y: " << cellX << ", " << cellY << std::endl;

    int particlesPerSide = (int)(sqrt(numParticles));
    // std::cout << "TRACER particlesPerSide: " << particlesPerSide << std::endl;

    real deltaX = grid_->dx() / (particlesPerSide);
    real deltaY = grid_->dy() / (particlesPerSide);
    // std::cout << "TRACER deltaX,Y: " << deltaX << ", " << deltaY << std::endl;

    for (int i = 1; i <= particlesPerSide; ++i)
    {
        for (int j = 1; j <= particlesPerSide; ++j)
        {
            Particle p(cellX + deltaX / 2 + (i - 1) * deltaX, cellY + deltaY / 2 + (j - 1) * deltaY, type);
            // std::cout << "TRACER pX,Y: " << cellX + deltaX / 2 + (i - 1) * deltaX << ", " << cellY + deltaY / 2 + (j - 1) * deltaY << std::endl;
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
    PROG("advancing particles");

    for (unsigned int i = 0; i < particles_.size(); ++i)
    {
        Particle *p = &particles_[i];
        real u = interpolateU(p->x(), p->y());
        real v = interpolateV(p->x(), p->y());

        // std::cout << "TRACER oldPosX,Y: " << p->x() << ", " << p->y() << std::endl;
        p->setX(p->x() + dt * u);
        p->setY(p->y() + dt * v);
        // std::cout << "TRACER newPosX,Y: " << p->x() << ", " << p->y() << std::endl;

        int newCellX = p->getCellX(grid_->dx());
        int newCellY = p->getCellY(grid_->dy());
        // std::cout << "TRACER cellX,Y: " << newCellX << ", " << newCellY << std::endl;

        // if the particle moved into an obstacle cell or outside the domain just delete it
        bool isOutsideDomain = newCellX < 1 || newCellX > grid_->imax() || newCellY < 1 || newCellY > grid_->jmax();
        bool isObstacle = grid_->isObstacle(newCellX, newCellY);
        // std::cout << "TRACER outside: " << isOutsideDomain << ", " << isObstacle << std::endl;

        if (isOutsideDomain || isObstacle)
        {
            particles_.erase(particles_.begin() + i);
            // the vector has now shrunk so we need to check the element which moved into the empty place
            i--;
        }
    }
}

real ParticleTracer::interpolateU(real x, real y)
{
    // see section 4.2.1
    // std::cout << "interpolateU x,y: " << x << ", " << y << std::endl;
    int i = (int)(x / grid_->dx()) + 1;
    int j = (int)((y + 0.5 * grid_->dy()) / grid_->dy()) + 1;
    // std::cout << "interpolateU i,j: " << i << ", " << j << std::endl;

    real x1, x2, y1, y2;
    x1 = (i - 1) * grid_->dx();
    x2 = i * grid_->dx();
    y1 = ((j - 1) - 0.5) * grid_->dy();
    y2 = (j - 0.5) * grid_->dy();
    // std::cout << "interpolateU x1,x2: " << x1 << ", " << x2 << std::endl;
    // std::cout << "interpolateU y1,y2: " << y1 << ", " << y2 << std::endl;

    // TODO: use u accessor function to incooperate obstacles
    real u1, u2, u3, u4;
    u1 = grid_->u()(i - 1 , j - 1 );
    u2 = grid_->u()(i     , j - 1 );
    u3 = grid_->u()(i - 1 , j     );
    u4 = grid_->u()(i     , j     );
    // std::cout << "interpolateU u1,u2: " << u1 << ", " << u2 << std::endl;
    // std::cout << "interpolateU u3,u4: " << u3 << ", " << u4 << std::endl;

    real u = (1 / (grid_->dx() * grid_->dy())) * (
                 (x2 - x) * (y2 - y) * u1 +
                 (x  - x1) * (y2 - y) * u2 +
                 (x2 - x) * (y  - y1) * u3 +
                 (x  - x1) * (y  - y1) * u4
             );

    return u;
}

real ParticleTracer::interpolateV(real x, real y)
{
    // see section 4.2.1

    // std::cout << "interpolateV x,y: " << x << ", " << y << std::endl;
    int i = (int)((x + 0.5 * grid_->dx()) / grid_->dx()) + 1;
    int j = (int)(y  / grid_->dy()) + 1;
    // std::cout << "interpolateV i,j: " << i << ", " << j << std::endl;

    real x1, x2, y1, y2;
    x1 = ((i - 1) - 0.5) * grid_->dx();
    x2 = (i - 0.5) * grid_->dx();
    y1 = (j - 1) * grid_->dy();
    y2 = j * grid_->dy();
    // std::cout << "interpolateV x1,x2: " << x1 << ", " << x2 << std::endl;
    // std::cout << "interpolateV y1,y2: " << y1 << ", " << y2 << std::endl;

    // TODO: use u accessor function to incooperate obstacles
    real v1, v2, v3, v4;
    v1 = grid_->v()(i - 1 , j - 1 );
    v2 = grid_->v()(i     , j - 1 );
    v3 = grid_->v()(i - 1 , j     );
    v4 = grid_->v()(i     , j     );

    real v = (1 / (grid_->dx() * grid_->dy())) * (
                 (x2 - x) * (y2 - y) * v1 +
                 (x  - x1) * (y2 - y) * v2 +
                 (x2 - x) * (y  - y1) * v3 +
                 (x  - x1) * (y  - y1) * v4
             );

    return v;
}


