#ifndef PARTICLE_HH
#define PARTICLE_HH

#include <iostream>
#include "Types.hh"


//*******************************************************************************************************************
/*    Particl class for 2 dimensions
*
*
*/
//*******************************************************************************************************************


class Particle
{

public:

    Particle(real x, real y);

    real &x()
    {
        return x_;
    }
    real &y()
    {
        return y_;
    }
    real &u()
    {
        return u_;
    }
    real &v()
    {
        return v_;
    }

    const real &x() const
    {
        return x_;
    }
    const real &y() const
    {
        return y_;
    }
    const real &u() const
    {
        return u_;
    }
    const real &v() const
    {
        return v_;
    }
    const int &type() const
    {
        return type_;
    }
    //inline real & operator () ( int i ,int j );
    //const  real & operator () ( int i ,int j ) const;

    int getCellX(real dx);
    int getCellY(real dy);

    void setX(real x);
    void setY(real y);

private:

    real x_ ;   // "x" position of particle
    real y_ ;   // "y" position of particle
    real u_ ;   // "u" velocity of particle
    real v_ ;   // "v" velocity of particle

    int type_;

};

inline void Particle::setX(real x) {
    x_ = x;
}
inline void Particle::setY(real y) {
    y_ = y;
}


#endif //PARTICLE_HHH
