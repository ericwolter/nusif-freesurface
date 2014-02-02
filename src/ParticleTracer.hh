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
	void fillCell(int i, int j, int numParticles, int type);
    void addRectangle(real x1, real y1, real x2, real y2, int type);
    void addCircle(real x, real y, real r, int type);

    void advanceParticles(real const dt);
    void particle_boundary(int i , int j ,real x,real y , real u , real v , const real dt);

    void print();

private:
    

    std::vector<Particle> particles_;
    StaggeredGrid *grid_;

    // Caution!!! These are only public for test purposes!!!
public:
    real interpolateU(real x, real y);
    real interpolateV(real x, real y);
};

#endif //PARTICLE_TRACER_HH
