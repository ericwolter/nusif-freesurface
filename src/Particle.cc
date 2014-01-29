#include "Particle.hh"

Particle::Particle(real xx, real yy, int type)
{
    x_ = xx;
    y_ = yy;
    type_ = type;
}

int Particle::getCellX(real dx)
{
    return (int)(x_ / dx) + 1;
}
int Particle::getCellY(real dy)
{
    return (int)(y_ / dy) + 1;
}
