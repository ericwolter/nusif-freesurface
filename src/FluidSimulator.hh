#ifndef __FLUID_SIMULATOR_H__
#define __FLUID_SIMULATOR_H__


#include "FileReader.hh"
#include "StaggeredGrid.hh"
#include "SORSolver.hh"

#include "VTKWriter.hh"

#include <cmath>


class FluidSimulator
{
  public:
      FluidSimulator( const FileReader & conf );

      /// Simulates a given time-length
      void simulate             ( real duration              );
      void simulateTimeStepCount( unsigned int nrOfTimeSteps );


      // Getter functions for the internally stored StaggeredGrid
            StaggeredGrid & grid() { return grid_; };
      const StaggeredGrid & grid() const { return grid_; };

            // Getter functions for the internally stored SORSolver
            SORSolver & solv() { return solver_; };
      const SORSolver & solv() const { return solver_; };

      
      void backstep_test();
      void normalization();
      void dcavity_test();
      void testFG();

  private:
      // helper functions (derivatives)
      real dxuu(int i, int j), dyuv(int i, int j), ddxu(int i, int j), ddyu(int i, int j);
      real dyvv(int i, int j), dxuv(int i, int j), ddxv(int i, int j), ddyv(int i, int j);

      void computeFG();
      void composeRHS();
      void updateVelocities();
      void determineNextDT( real const & limit );
      void refreshBoundaries();

      // needed values
      real safetyfac_, gamma_, Re_, gx_, gy_, dt_, vel_N, vel_S, vel_E, vel_W, uInit_, vInit_, pInit_;
      real rectX_, rectXX_, rectY_, rectYY_, circX_, circY_, circR_;
      BCTYPE cond_N, cond_S, cond_E, cond_W;
      int timeStepNr, normfreq, outPutInt, imax, jmax;

  protected:
      StaggeredGrid grid_;   //< grid
      SORSolver solver_;     //< solver

};



#endif
