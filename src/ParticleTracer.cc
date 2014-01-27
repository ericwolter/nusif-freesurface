#include "math.h"
#include <iostream>

#include "ParticleTracer.hh"


ParticleTracer::ParticleTracer ( const FileReader &conf  )
    : grid_(conf)
{
}

ParticleTracer::ParticleTracer ( const StaggeredGrid &grid )
    : grid_(grid)
{
}

void ParticleTracer::markCells()
{
    for (int i = 1; i <= grid_.imax(); ++i)
    {
        for (int j = 1; j <= grid_.jmax(); ++j)
        {
            if (!grid_.isFluid(i, j)) continue;

            grid_.setCellToEmpty(i, j);
        }
    }

    for (std::vector<Particle>::iterator p = particles_.begin() ; p != particles_.end(); ++p)
    {
        int i = p->getCellX(grid_.dx());
        int j = p->getCellY(grid_.dy());

        grid_.setCellToFluid(i, j);
    }
}

void ParticleTracer::addRectangle(int x1, int y1, int x2, int y2)
{
    std::cout << "Adding particle rectangle: " << "(" << x1 << "|" << y1 << ")" << ", " << "(" << x2 << "|" << y2 << ")" << std::endl;
    int minX = std::min(x1, x2);
    int maxX = std::max(x1, x2);
    int minY = std::min(y1, y2);
    int maxY = std::max(y1, y2);

    for (int i = minX; i < maxX; ++i)
    {
        for (int j = minY; j < maxY; ++j)
        {
            // std::cout << "TRACER fillCell: " << i + 1 << ", " << j + 1 << std::endl;
            this->fillCell(i + 1, j + 1, 9);
        }
    }
}

void ParticleTracer::addCircle(int x, int y, int r)
{
    std::cout << "Adding particle circle: " << "(" << x << "|" << y << "|" << r << ")" << std::endl;

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
    real cellX = (x - 1) * grid_.dx();
    real cellY = (y - 1) * grid_.dy();
    // std::cout << "TRACER cellX,Y: " << cellX << ", " << cellY << std::endl;

    int particlesPerSide = (int)(sqrt(numParticles));
    // std::cout << "TRACER particlesPerSide: " << particlesPerSide << std::endl;

    real deltaX = grid_.dx() / (particlesPerSide + 1);
    real deltaY = grid_.dy() / (particlesPerSide + 1);
    // std::cout << "TRACER deltaX,Y: " << deltaX << ", " << deltaY << std::endl;

    for (int i = 1; i <= particlesPerSide; ++i)
    {
        for (int j = 1; j <= particlesPerSide; ++j)
        {
            Particle p(cellX + i * deltaX, cellY + j * deltaY);
            // std::cout << "TRACER pX,Y: " << (cellX + i * deltaX) << ", " << (cellY + j * deltaY) << std::endl;
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
    // see section 4.2.1

    // std::cout << "interpolateU x,y: " << x << ", " << y << std::endl;
    int i = (int) (x / grid_.dx()) + 1;
    int j = (int) (( y + 0.5 * grid_.dy() ) / grid_.dy()) + 1;
    // std::cout << "interpolateU i,j: " << i << ", " << j << std::endl;

    real x1, x2, y1, y2;
    x1 = (i - 1) * grid_.dx();
    x2 = i * grid_.dx();
    y1 = ((j - 1) - 0.5) * grid_.dy();
    y2 = (j - 0.5) * grid_.dy();
    // std::cout << "interpolateU x1,x2: " << x1 << ", " << x2 << std::endl;
    // std::cout << "interpolateU y1,y2: " << y1 << ", " << y2 << std::endl;

    // TODO: use u accessor function to incooperate obstacles
    real u1, u2, u3, u4;
    u1 = grid_.u()( i - 1 , j - 1 );
    u2 = grid_.u()( i     , j - 1 );
    u3 = grid_.u()( i - 1 , j     );
    u4 = grid_.u()( i     , j     );

    real u = (1 / (grid_.dx() * grid_.dy())) * (
                 ( x2 - x  ) * ( y2 - y  ) * u1 +
                 ( x  - x1 ) * ( y2 - y  ) * u2 +
                 ( x2 - x  ) * ( y  - y1 ) * u3 +
                 ( x  - x1 ) * ( y  - y1 ) * u4
             );

    return u;
}

real ParticleTracer::interpolateV(real x, real y)
{
    // see section 4.2.1

    // std::cout << "interpolateV x,y: " << x << ", " << y << std::endl;
    int i = (int) (( x + 0.5 * grid_.dx() ) / grid_.dx()) + 1;
    int j = (int) ( y  / grid_.dy()) + 1;
    // std::cout << "interpolateV i,j: " << i << ", " << j << std::endl;

    real x1, x2, y1, y2;
    x1 = ((i - 1) - 0.5) * grid_.dx();
    x2 = (i - 0.5) * grid_.dx();
    y1 = (j - 1) * grid_.dy();
    y2 = j * grid_.dy();
    // std::cout << "interpolateV x1,x2: " << x1 << ", " << x2 << std::endl;
    // std::cout << "interpolateV y1,y2: " << y1 << ", " << y2 << std::endl;

    // TODO: use u accessor function to incooperate obstacles
    real v1, v2, v3, v4;
    v1 = grid_.v()( i - 1 , j - 1 );
    v2 = grid_.v()( i     , j - 1 );
    v3 = grid_.v()( i - 1 , j     );
    v4 = grid_.v()( i     , j     );

    real v = (1 / (grid_.dx() * grid_.dy())) * (
                 ( x2 - x  ) * ( y2 - y  ) * v1 +
                 ( x  - x1 ) * ( y2 - y  ) * v2 +
                 ( x2 - x  ) * ( y  - y1 ) * v3 +
                 ( x  - x1 ) * ( y  - y1 ) * v4
             );

    return v;
}


