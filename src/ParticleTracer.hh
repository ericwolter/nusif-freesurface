#ifndef PARTICLE_TRACER_HH
#define PARTICLE_TRACER_HH

#include <vector>

#include "Particle.hh"
#include "StaggeredGrid.hh"

class ParticleTracer
{
public:

    ParticleTracer ( StaggeredGrid &grid );

    const std::vector<Particle> &particles() const
    {
        return particles_;
    }

    void markCells();
    void addRectangle(int x1, int y1, int x2, int y2);
    void addCircle(int x, int y, int r);

    void print();
private:
    void fillCell(int i, int j, int numParticles);

    std::vector<Particle> particles_;
    StaggeredGrid grid_;
};

#endif //PARTICLE_TRACER_HH
