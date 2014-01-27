#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"
#include "ParticleTracer.hh"

#include "VTKWriter.hh"

#include <cmath>


class FluidSimulator
{
public:
    FluidSimulator( const FileReader &conf );

    /// Simulates a given time-length
    void simulate             ( real duration              );
    void simulateTimeStepCount( unsigned int nrOfTimeSteps );


    // Getter functions for the internally stored StaggeredGrid
    StaggeredGrid &grid()
    {
        return grid_;
    };
    const StaggeredGrid &grid() const
    {
        return grid_;
    };

    // Getter functions for the internally stored SORSolver
    SORSolver &solv()
    {
        return solver_;
    };
    const SORSolver &solv() const
    {
        return solver_;
    };


    void backstep_test();
    void normalization();
    void dcavity_test();
    void testFG();

    // Compute surface boundary
    void set_UVP_surface(const int &dt);
    void set_UVP_surface(int i, int j , const int &dt);

    void one_empty_neighbour   (int i , int j , const int &dt ) ;
    void two_empty_neighbour   (int i , int j , const int &dt ) ;
    void three_empty_neighbour (int i , int j , const int &dt ) ;
    void four_empty_neighbour  (int i , int j , const int &dt ) ;
private:
    // helper functions (derivatives)
    real dxuu(int i, int j), dyuv(int i, int j), ddxu(int i, int j), ddyu(int i, int j);
    real dyvv(int i, int j), dxuv(int i, int j), ddxv(int i, int j), ddyv(int i, int j);

    void computeFG();
    void composeRHS();
    void updateVelocities();
    void determineNextDT( real const &limit );
    void refreshBoundaries();

    //void set_UVP_surface(const int& dt);
    // needed values
    real safetyfac_, gamma_, Re_, gx_, gy_, dt_, vel_N, vel_S, vel_E, vel_W, uInit_, vInit_, pInit_;
    real rectX_, rectXX_, rectY_, rectYY_, circX_, circY_, circR_;
    BCTYPE cond_N, cond_S, cond_E, cond_W;
    unsigned int timeStepNr, normfreq, outPutInt;
    int imax, jmax;
    int rectX1_particle_ , rectX2_particle_ , rectY1_particle_ , rectY2_particle_ ;

protected:
    StaggeredGrid grid_;   //< grid
    SORSolver solver_;     //< solver
    ParticleTracer particle_tracer_ ;
};



#endif
