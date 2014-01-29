#ifndef PARTICLE_TRACER_HH
#define PARTICLE_TRACER_HH

#include <vector>

#include "Particle.hh"
#include "StaggeredGrid.hh"

class ParticleTracer
{
public:

    ParticleTracer();
    ParticleTracer(StaggeredGrid *grid);

    const std::vector<Particle> &particles() const
    {
        return particles_;
    }

    void markCells();
    void addRectangle(int x1, int y1, int x2, int y2, int type);
    void addCircle(int x, int y, int r, int type);

    void advanceParticles(real const dt);
    void print();

private:
    void fillCell(int i, int j, int numParticles, int type);

    std::vector<Particle> particles_;
    StaggeredGrid *grid_;

    // Caution!!! These are only public for test purposes!!!
public:
    real interpolateU(real x, real y);
    real interpolateV(real x, real y);
};

#endif //PARTICLE_TRACER_HH
