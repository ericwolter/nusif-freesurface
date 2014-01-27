#include "Particle.hh"

Particle::Particle(real xx, real yy)
{
    x_ = xx;
    y_ = yy;
}

int Particle::getCellX(real dx)
{
    return (int)(x_ / dx);
}
int Particle::getCellY(real dy)
{
    return (int)(y_ / dy);
}
