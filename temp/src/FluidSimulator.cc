
#include "FluidSimulator.hh"
#include <iostream>



// constructor for FluidSimulator
FluidSimulator::FluidSimulator( const FileReader &conf )
    : grid_(conf), solver_(conf)
{
    PROG("construct FluidSimulator with a file");
    gamma_ = conf.getRealParameter("gamma");
    CHECK_MSG( (gamma_ >= 0 && gamma_ <= 1), "wrong input for gamma: " << gamma_);
    Re_ = conf.getRealParameter("Re");
    CHECK_MSG( (Re_ > 0), "wrong input for Re: " << Re_);
    gx_ = conf.getRealParameter("GX");
    gy_ = conf.getRealParameter("GY");
    dt_ = conf.getRealParameter("dt");
    CHECK_MSG( (dt_ > 0), "wrong input for dt: " << dt_);

    vel_N = conf.getRealParameter("boundary_velocity_N");
    vel_S = conf.getRealParameter("boundary_velocity_S");
    vel_E = conf.getRealParameter("boundary_velocity_E");
    vel_W = conf.getRealParameter("boundary_velocity_W");

    std::string no = conf.getStringParameter("boundary_condition_N");
    std::string so = conf.getStringParameter("boundary_condition_S");
    std::string we = conf.getStringParameter("boundary_condition_W");
    std::string es = conf.getStringParameter("boundary_condition_E");

    if ( no == "inflow" )
    {
        cond_N = INFLOW;
    }
    else if ( no == "outflow" )
    {
        ASSERT_MSG(vel_N == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_N = OUTFLOW;
    }
    else
    {
        cond_N = NOSLIP;
    }

    if ( so == "inflow" )
    {
        cond_S = INFLOW;
    }
    else if ( so == "outflow" )
    {
        ASSERT_MSG(vel_S == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_S = OUTFLOW;
    }
    else
    {
        cond_S = NOSLIP;
    }

    if ( we == "inflow" )
    {
        cond_W = INFLOW;
    }
    else if ( we == "outflow" )
    {
        ASSERT_MSG(vel_W == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_W = OUTFLOW;
    }
    else
    {
        cond_W = NOSLIP;
    }

    if ( es == "inflow" )
    {
        cond_E = INFLOW;
    }
    else if ( es == "outflow" )
    {
        ASSERT_MSG(vel_E == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_E = OUTFLOW;
    }
    else
    {
        cond_E = NOSLIP;
    }

    safetyfac_ = conf.getRealParameter("safetyfactor");
    timeStepNr = (unsigned int)conf.getIntParameter("timesteps");
    CHECK_MSG( (timeStepNr >= 0), "wrong input for timesteps: " << timeStepNr);
    outPutInt = (unsigned int)conf.getIntParameter("outputinterval");
    CHECK_MSG( (outPutInt > 0), "wrong input for outputinterval: " << outPutInt);
    normfreq = (unsigned int)conf.getIntParameter("normalizationfrequency");
    CHECK_MSG( (normfreq > 0), "wrong input for normalizationfrequency: " << normfreq);
    uInit_ = conf.getRealParameter("U_INIT");
    vInit_ = conf.getRealParameter("V_INIT");
    pInit_ = conf.getRealParameter("P_INIT");
    imax = conf.getIntParameter("imax");
    CHECK_MSG( (imax > 0), "wrong input for imax: " << imax);
    jmax = conf.getIntParameter("jmax");
    CHECK_MSG( (jmax > 0), "wrong input for jmax: " << jmax);

    // obstacles:
    // rectangle
    rectX_ = conf.getRealParameter("RectangleX1");
    CHECK_MSG( (rectX_ >= 0), "wrong input for RectangleX1: " << rectX_);
    rectXX_ = conf.getRealParameter("RectangleX2");
    CHECK_MSG( (rectXX_ >= 0), "wrong input for RectangleX2: " << rectXX_);
    rectY_ = conf.getRealParameter("RectangleY1");
    CHECK_MSG( (rectY_ >= 0), "wrong input for RectangleY1: " << rectY_);
    rectYY_ = conf.getRealParameter("RectangleY2");
    CHECK_MSG( (rectYY_ >= 0), "wrong input for RectangleY2: " << rectYY_);
    // circle
    circX_ = conf.getRealParameter("CircleX");
    CHECK_MSG( (circX_ >= 0), "wrong input for CircleX: " << circX_);
    circY_ = conf.getRealParameter("CircleY");
    CHECK_MSG( (circY_ >= 0), "wrong input for CircleY: " << circY_);
    circR_ = conf.getRealParameter("CircleR");
    CHECK_MSG( (circR_ >= 0), "wrong input for CircleR: " << circR_);

}

void FluidSimulator::composeRHS()
{
    // formula 3.38
    for ( int i = 1; i <= imax; i++ )
    {
        for ( int j = 1; j <= jmax; j++ )
            grid_.rhs()(i - 1, j - 1) = ( (grid_.f()(i, j - 1) - grid_.f()(i - 1, j - 1)) / grid_.dx() + (grid_.g()(i - 1, j) - grid_.g()(i - 1, j - 1)) / grid_.dy() ) / dt_;
    }
}

void FluidSimulator::updateVelocities()
{
    real dtdx = dt_ / grid_.dx();
    real dtdy = dt_ / grid_.dy();

    for ( int i = 1; i <= imax; i++ )
    {
        for ( int j = 1; j <= jmax; j++ )
        {

            if ( grid_.isFluid(i, j) )  // update obly for fluid cells
            {

                if ( i != imax && grid_.isFluid(i + 1, j) ) // formula 3.34
                    grid_.u()(i, j) = grid_.f()(i, j - 1) - ( grid_.p(i + 1, j, WEST) - grid_.p()(i, j) ) * dtdx;

                if ( j != jmax && grid_.isFluid(i, j + 1) ) // formula 3.35
                    grid_.v()(i, j) = grid_.g()(i - 1, j) - ( grid_.p(i, j + 1, SOUTH) - grid_.p()(i, j) ) * dtdy;

            }
        }
    }
}

void FluidSimulator::determineNextDT( real const &limit )
{
    // first stability condition
    real one = Re_ / ( 2 / pow(grid_.dx(), 2) + 2 / pow(grid_.dy(), 2) );
    bool first = dt_ < one;
    // CFL / 3.49
    real two = grid_.dx() / grid_.u().maxE();
    bool second = dt_ < two;
    real three = grid_.dy() / grid_.v().maxE();
    bool third = dt_ < three;
    // control step size
    if ( !(first && second && third) )   // 3.50
    {
        real minimum = std::min(one, two);
        minimum = std::min(minimum, three);

        if ( limit > 0 && limit <= 1 )   // check safety factor
        {
            dt_ = limit * minimum;
        }
        else
        {
            dt_ = minimum * 0.5;
        }
        PROG("dt is now: " << dt_);
    }

}

void FluidSimulator::refreshBoundaries()
{
    // conditions for north border
    if ( cond_N == NOSLIP )   // no-slip
    {
        for ( int i = 1; i <= imax; i++ )
        {
            grid_.v()(i, jmax) = 0;
            grid_.u()(i, jmax + 1) = 2 * vel_N - grid_.u()(i, jmax);
        }

    }
    else if ( cond_N == INFLOW )     // inflow
    {
        for ( int i = 1; i <= imax; i++ )
        {
            if ( grid_.isFluid(i, jmax) )
            {
                grid_.v()(i, jmax) = vel_N;
                grid_.u()(i, jmax + 1) = - grid_.u()(i, jmax);
            }
            else
            {
                grid_.v()(i, jmax) = 0;
                grid_.u()(i, jmax + 1) = 2 * vel_N - grid_.u()(i, jmax);
            }
        }

    }
    else     // outflow
    {
        for ( int i = 1; i <= imax; i++ )
        {
            grid_.v()(i, jmax) = grid_.v()(i, jmax - 1);
            grid_.u()(i, jmax + 1) = grid_.u()(i, jmax);
        }
    }

    // conditions for south border
    if ( cond_S == NOSLIP )   // no-slip
    {
        for ( int i = 1; i <= imax; i++ )
        {
            grid_.v()(i, 0) = 0;
            grid_.u()(i, 0) = 2 * vel_S - grid_.u()(i, 1);
        }

    }
    else if ( cond_S == INFLOW )     // inflow
    {
        for ( int i = 1; i <= imax; i++ )
        {
            if ( grid_.isFluid(i, 1) )
            {
                grid_.v()(i, 0) = vel_S;
                grid_.u()(i, 0) = - grid_.u()(i, 1);
            }
            else
            {
                grid_.v()(i, 0) = 0;
                grid_.u()(i, 0) = 2 * vel_S - grid_.u()(i, 1);
            }
        }

    }
    else     // outflow
    {
        for ( int i = 1; i <= imax; i++ )
        {
            grid_.v()(i, 0) = grid_.v()(i, 1);
            grid_.u()(i, 0) = grid_.u()(i, 1);
        }
    }

    // conditions for west border
    if ( cond_W == NOSLIP )   // no-slip
    {
        for ( int j = 1; j <= jmax; j++ )
        {
            grid_.u()(0, j) = 0;
            grid_.v()(0, j) = 2 * vel_W - grid_.v()(1, j);
        }

    }
    else if ( cond_W == INFLOW )     // inflow
    {
        for ( int j = 1; j <= jmax; j++ )
        {
            if ( grid_.isFluid(1, j) )
            {
                grid_.u()(0, j) = vel_W;
                grid_.v()(0, j) = - grid_.v()(1, j);
            }
            else
            {
                grid_.u()(0, j) = 0;
                grid_.v()(0, j) = 2 * vel_W - grid_.v()(1, j);
            }

        }

    }
    else     // outflow
    {
        for ( int j = 1; j <= jmax; j++ )
        {
            grid_.u()(0, j) = grid_.u()(1, j);
            grid_.v()(0, j) = grid_.v()(1, j);
        }
    }

    // conditions for east border
    if ( cond_E == NOSLIP )   // no-slip
    {
        for ( int j = 1; j <= jmax; j++ )
        {
            grid_.u()(imax, j) = 0;
            grid_.v()(imax + 1, j) = 2 * vel_E - grid_.v()(imax, j);
        }

    }
    else if ( cond_E == INFLOW )     // inflow
    {
        for ( int j = 1; j <= jmax; j++ )
        {
            if ( grid_.isFluid(imax, j) )
            {
                grid_.u()(imax, j) = vel_E;
                grid_.v()(imax + 1, j) = - grid_.v()(imax, j);
            }
            else
            {
                grid_.u()(imax, j) = 0;
                grid_.v()(imax + 1, j) = 2 * vel_E - grid_.v()(imax, j);
            }

        }

    }
    else     // outflow
    {
        for ( int j = 1; j <= jmax; j++ )
        {
            grid_.u()(imax, j) = grid_.u()(imax - 1, j);
            grid_.v()(imax + 1, j) = grid_.v()(imax, j);
        }
    }

}

void FluidSimulator::simulate( real duration )
{
    VTKWriter vtkWriter ( "lidDrivenCavity" );
    real t = 0;
    unsigned int n = 0;

    PROG("initialize u, v, p, rhs");
    int half = (int) ( rectYY_ / grid_.dy() );
    for ( int i = 0; i < grid_.u().getSize(0); ++i )
    {
        for ( int j = 0; j < grid_.u().getSize(1); ++j )
        {
            if ( j > half )
            {
                grid_.u()(i, j) = 1;
            }
            else
            {
                grid_.u()(i, j) = 0;
            }
        }
    }

    grid_.v().fill(vInit_);
    grid_.p().fill(pInit_);
    grid_.rhs().fill(0);
    PROG("set inner obstacles");
    grid_.createRectangle(rectX_, rectY_, rectXX_, rectYY_);
    grid_.createCircle(circX_, circY_, circR_);

    refreshBoundaries();

    while (t <= duration )
    {
        // without 0th step:
        // if ( n%outPutInt == 0 && n != 0 )
        if ( n % outPutInt == 0 )
            vtkWriter.write(grid_, NULL);
        PROG(n << "'th timestep: determine next dt");
        determineNextDT( safetyfac_ );
        PROG(n << "'th timestep: refresh boundaries");
        refreshBoundaries();
        computeFG();
        PROG(n << "'th timestep: compute the right-hand side of the pressure equation");
        composeRHS();
        PROG(n << "'th timestep: solve pressure equation");
        solv().solve( grid_ );
        PROG(n << "'th timestep: update u and v; " );
        updateVelocities();
        if ( n % normfreq == 0 )
            normalization();
        t += dt_;
        n++;
    }
}

void FluidSimulator::simulateTimeStepCount( unsigned int nrOfTimeSteps )
{
    VTKWriter vtkWriter ("lidDrivenCavity");
    unsigned int n = 0;

    PROG("initialize u, v, p, rhs");
    int half = (int) ( rectYY_ / grid_.dy() );
    for ( int i = 0; i < grid_.u().getSize(0); ++i )
    {
        for ( int j = 0; j < grid_.u().getSize(1); ++j )
        {
            if ( j > half )
            {
                grid_.u()(i, j) = 1;
            }
            else
            {
                grid_.u()(i, j) = 0;
            }
        }
    }

    grid_.v().fill(vInit_);
    grid_.p().fill(pInit_);
    grid_.rhs().fill(0);
    PROG("set inner obstacles");
    grid_.createRectangle(rectX_, rectY_, rectXX_, rectYY_);
    grid_.createCircle(circX_, circY_, circR_);
    grid_.createPng("test.png");

    refreshBoundaries();

    while (n <= nrOfTimeSteps )
    {
        // without 0th step:
        // if ( n%outPutInt == 0 && n != 0 )
        if ( n % outPutInt == 0 )
            vtkWriter.write(grid_, NULL);
        PROG(n << "'th timestep: determine next dt");
        determineNextDT( safetyfac_ );
        PROG(n << "'th timestep: refresh boundaries");
        refreshBoundaries();
        computeFG();
        PROG(n << "'th timestep: compute the right-hand side of the pressure equation");
        composeRHS();
        PROG(n << "'th timestep: solve pressure equation");
        solv().solve( grid_ );
        PROG(n << "'th timestep: update u and v");
        updateVelocities();
        if ( n % normfreq == 0 )
            normalization();
        n++;
    }

}

void FluidSimulator::normalization()
{
    real psum = 0;

    for ( int i = 1; i <= imax; i++ )
    {
        for ( int j = 1; j <= jmax; j++ )
        {

            if ( grid_.isFluid(i, j) )
                psum += grid_.p()(i, j);
        }
    }

    real pMean = psum / grid_.getNumFluid();
    // substract the mean
    for ( int i = 1; i <= imax; i++ )
    {
        for ( int j = 1; j <= jmax; j++ )
        {

            if ( grid_.isFluid(i, j) )
                grid_.p()(i, j) -= pMean;
        }
    }

}

// help functions for derivatives
// d(u^2)/dx
real FluidSimulator::dxuu(int i, int j)
{
    ASSERT_MSG( (i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG( (i + 1 < grid_.u().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG( (j >= 0), "wrong input for j: " << j);
    real res = (grid_.u()(i, j) + grid_.u(i + 1, j, WEST)) * (grid_.u()(i, j) + grid_.u(i + 1, j, WEST)) / (4 * grid_.dx());
    res -= (grid_.u(i - 1, j, EAST) + grid_.u()(i, j)) * (grid_.u(i - 1, j, EAST) + grid_.u()(i, j)) / (4 * grid_.dx());
    res += gamma_ * ( fabs(grid_.u()(i, j) + grid_.u(i + 1, j, WEST)) * (grid_.u()(i, j) - grid_.u(i + 1, j, WEST)) * 0.25 ) / grid_.dx();
    res -= gamma_ * ( fabs(grid_.u(i - 1, j, EAST) + grid_.u()(i, j)) * (grid_.u(i - 1, j, EAST) - grid_.u()(i, j)) * 0.25 ) / grid_.dx();

    return res;
}

// d(v^2)/dy
real FluidSimulator::dyvv(int i, int j)
{
    ASSERT_MSG( (i >= 0), "wrong input for i: " << i);
    ASSERT_MSG( (j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG( (j + 1 < grid_.v().getSize(1)), "wrong input for j: " << j);
    real res = (grid_.v()(i, j) + grid_.v(i, j + 1, SOUTH)) * (grid_.v()(i, j) + grid_.v(i, j + 1, SOUTH)) / (4 * grid_.dy());
    res -= (grid_.v(i, j - 1, NORTH) + grid_.v()(i, j)) * (grid_.v(i, j - 1, NORTH) + grid_.v()(i, j)) / (4 * grid_.dy());
    res += gamma_ * ( fabs(grid_.v()(i, j) + grid_.v(i, j + 1, SOUTH)) * (grid_.v()(i, j) - grid_.v(i, j + 1, SOUTH)) * 0.25 ) / grid_.dy();
    res -= gamma_ * ( fabs(grid_.v(i, j - 1, NORTH) + grid_.v()(i, j)) * (grid_.v(i, j - 1, NORTH) - grid_.v()(i, j)) * 0.25 ) / grid_.dy();

    return res;
}
////////////////////////////////////////////////////////////////////////////////
// d(uv)/dy
real FluidSimulator::dyuv(int i, int j)
{
    ASSERT_MSG( (i >= 0), "wrong input for i: " << i);
    ASSERT_MSG( (i + 1 < grid_.v().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG( (j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG( (j + 1 < grid_.u().getSize(1)), "wrong input for j: " << j);
    real diag = 0.0;
    if ( grid_.isFluid(i + 1, j - 1) )
    {
        diag = grid_.v()(i + 1, j - 1);
    }
    else
    {
        diag = 0.5 * ( grid_.v(i + 1, j, WEST) + grid_.v(i, j - 1, NORTH) );
    }
    real res = ( (grid_.v()(i, j) + grid_.v(i + 1, j, WEST)) * (grid_.u()(i, j) + grid_.u(i, j + 1, SOUTH)) * 0.25 ) / grid_.dy();
    res -= ( (grid_.v(i, j - 1, NORTH) + diag) * (grid_.u(i, j - 1, NORTH) + grid_.u()(i, j)) * 0.25 ) / grid_.dy();
    res += gamma_ * ( fabs(grid_.v()(i, j) + grid_.v(i + 1, j, WEST)) * (grid_.u()(i, j) - grid_.u(i, j + 1, SOUTH)) * 0.25 ) / grid_.dy();
    res -= gamma_ * ( fabs(grid_.v(i, j - 1, NORTH) + diag) * (grid_.u(i, j - 1, NORTH) - grid_.u()(i, j)) * 0.25 ) / grid_.dy();

    return res;
}

// d(uv)/dx
real FluidSimulator::dxuv(int i, int j)
{
    ASSERT_MSG( (j >= 0), "wrong input for j: " << j);
    ASSERT_MSG( (j + 1 < grid_.u().getSize(1)), "wrong input for j: " << j);
    ASSERT_MSG( (i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG( (i + 1 < grid_.v().getSize(0)), "wrong input for i: " << i);
    real diag = 0.0;
    if ( grid_.isFluid(i - 1, j + 1) )
    {
        diag = grid_.u()(i - 1, j + 1);
    }
    else
    {
        diag = 0.5 * ( grid_.u(i, j + 1, EAST) + grid_.u(i - 1, j, SOUTH) );
    }
    real res = ( (grid_.u()(i, j) + grid_.u(i, j + 1, SOUTH)) * (grid_.v()(i, j) + grid_.v(i + 1, j, WEST)) * 0.25 ) / grid_.dx();
    res -= ( (grid_.u(i - 1, j, EAST) + diag) * (grid_.v(i - 1, j, EAST) + grid_.v()(i, j)) * 0.25 ) / grid_.dx();
    res += gamma_ * ( fabs(grid_.u()(i, j) + grid_.u(i, j + 1, SOUTH)) * (grid_.v()(i, j) - grid_.v(i + 1, j, WEST)) * 0.25 ) / grid_.dx();
    res -= gamma_ * ( fabs(grid_.u(i - 1, j, EAST) + diag) * (grid_.v(i - 1, j, EAST) - grid_.v()(i, j)) * 0.25 ) / grid_.dx();

    return res;
}

// d^2u/dx^2
real FluidSimulator::ddxu(int i, int j)
{
    ASSERT_MSG( (i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG( (i + 1 < grid_.u().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG( (j >= 0), "wrong input for j: " << j);
    real res = ( grid_.u(i + 1, j, WEST) - 2 * grid_.u()(i, j) + grid_.u(i - 1, j, EAST) ) / (grid_.dx() * grid_.dx());

    return res;
}

// d^2v/dx^2
real FluidSimulator::ddxv(int i, int j)
{
    ASSERT_MSG( (i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG( (i + 1 < grid_.v().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG( (j >= 0), "wrong input for j: " << j);
    real res = ( grid_.v(i + 1, j, WEST) - 2 * grid_.v()(i, j) + grid_.v(i - 1, j, EAST) ) / (grid_.dx() * grid_.dx());

    return res;
}

// d^2u/dy^2
real FluidSimulator::ddyu(int i, int j)
{
    ASSERT_MSG( (j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG( (j + 1 < grid_.u().getSize(1)), "wrong input for j: " << j);
    ASSERT_MSG( (i >= 0), "wrong input for i: " << i);
    real res = ( grid_.u(i, j + 1, SOUTH) - 2 * grid_.u()(i, j) + grid_.u(i, j - 1, NORTH) ) / (grid_.dy() * grid_.dy());

    return res;
}

// d^2v/dy^2
real FluidSimulator::ddyv(int i, int j)
{
    ASSERT_MSG( (j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG( (j + 1 < grid_.v().getSize(1)), "wrong input for j: " << j);
    ASSERT_MSG( (i >= 0), "wrong input for i: " << i);
    real res = ( grid_.v(i, j + 1, SOUTH) - 2 * grid_.v()(i, j) + grid_.v(i, j - 1, NORTH) ) / (grid_.dy() * grid_.dy());

    return res;
}

void FluidSimulator::computeFG()
{
    grid_.f().fill(0);
    grid_.g().fill(0);
    PROG("set boundary values for F");
    // boundary values for F
    for ( int j = 1; j <= jmax; j++ )
    {
        grid_.f()(0, j - 1) = grid_.u()(0, j);
        grid_.f()(imax, j - 1) = grid_.u()(imax, j);
    }

    PROG("set boundary values for G");
    // boundary values for G
    for ( int i = 1; i <= imax; i++ )
    {
        grid_.g()(i - 1, 0) = grid_.v()(i, 0);
        grid_.g()(i - 1, jmax) = grid_.v()(i, jmax);
    }

    real dtRe = dt_ / Re_;
    real dtGx = dt_ * gx_;
    real dtGy = dt_ * gy_;

    PROG("compute F");
    // compute F
    for ( int i = 1; i < imax; i++ )
    {
        for ( int j = 1; j <= jmax; j++ )
        {

            if ( grid_.isFluid(i, j) && grid_.isFluid(i + 1, j) )
                grid_.f()(i, j - 1) = grid_.u()(i, j) + dtRe * (ddxu(i, j) + ddyu(i, j)) - dt_ * (dxuu(i, j) + dyuv(i, j)) + dtGx;

        }
    }

    PROG("compute G");
    // compute G
    for ( int i = 1; i <= imax; i++ )
    {
        for ( int j = 1; j < jmax; j++ )
        {

            if ( grid_.isFluid(i, j) && grid_.isFluid(i, j + 1) )
                grid_.g()(i - 1, j) = grid_.v()(i, j) + dtRe * (ddxv(i, j) + ddyv(i, j)) - dt_ * (dxuv(i, j) + dyvv(i, j)) + dtGy;

        }
    }

}


// helper functions for testing

void FluidSimulator::backstep_test()
{
    //simulate(dt_*timeStepNr);
    simulateTimeStepCount(timeStepNr);
}

void FluidSimulator::dcavity_test()
{
    //simulate(dt_*timeStepNr);
    simulateTimeStepCount(timeStepNr);
}

void FluidSimulator::testFG()
{
    computeFG();
    //grid_.rhs().print();
}

//*******************************************************************************************************************

void FluidSimulator::set_UVP_surface(const int& dt)
{

int number_of_empty_neighbour = 0 ;
    for ( int i = 1; i <= imax; i++ )
    {
        for ( int j = 1; j < jmax; j++ ) // !?!
        {

            if ( grid_.isFluid(i, j) )
            	{
            		// Count the number of empty neighbour cell
            		if ( grid_.isEmpty(i+1, j) )
            			++ number_of_empty_neighbour ;
            		if ( grid_.isEmpty(i-1, j) )
        			++ number_of_empty_neighbour ;
            		if ( grid_.isEmpty(i, j+1) )
            			++ number_of_empty_neighbour ;
            		if ( grid_.isEmpty(i, j-1) )
            			++ number_of_empty_neighbour ;

            	}
            switch (number_of_empty_neighbour )
            {
            		case(1):
            			one_empty_neighbour  ( i , j , dt ) ;
            		case(2):
            			two_empty_neighbour  ( i , j , dt ) ;
            		case(3):
            			three_empty_neighbour( i , j , dt ) ;
            		case(4):
            			four_empty_neighbour ( i , j , dt ) ;
            }
        }
        }
    }



//*******************************************************************************************************************

void FluidSimulator::one_empty_neighbour(int i , int j , const int& dt)
{
	// According to page 92 ,case 1 of the book
	// If cell (i,j) is Fluid, one of the other neighbour(North,East,West,South) is empty cell.
	// In any case of empty neighbour cell the formula is the same but we need to change for indices for mentioning
	// right neighbour cell.

	// EAST : only the cell (i+1,j) is empty (exactly like one which has been mentioned in the book)

	if ( grid_.isEmpty(i+1,j) )
	 {
		grid_.p()(i,j) =  ( 2 / (Re_*grid_.dx()) ) * ( grid_.u()( i , j) - grid_.u( i-1 , j , EAST) ) ;
		grid_.u()(i,j) =  grid_.u(i-1 , j , EAST) - ( grid_.dx()/ grid_.dy()) * ( grid_.v()( i , j ) - grid_.v( i , j-1 , NORTH) ) ;

		if (grid_.isEmpty(i+1,j-1))
			grid_.v()(i+1,j-1) = grid_.v(i , j-1 , NORTH ) - ( grid_.dx()/ grid_.dy()) * ( grid_.u()( i , j ) - grid_.u( i , j-1 , NORTH) ) ;
		else ///////////////////////////////////////////// PROBLEM ////////////////////////////////////////////////////////////////////////
			grid_.v()(i+1,j-1) = grid_.v(i+1 , j-2 , NORTH ) - ( grid_.dy()/ grid_.dx()) * ( grid_.u( i+1 , j-1 , NORTH ) - grid_.u( i , j-1 , NORTH) ) ;

	 }

	// WEST : only the cell (i-1,j) is empty (West of cell(i,j) )

	else if ( grid_.isEmpty(i-1,j) )
	 {
		grid_.p()(i,j) =  ( 2 / (Re_*grid_.dx()) ) * ( grid_.u( i+1 , j , WEST) - grid_.u()( i , j) ) ;
		grid_.u()(i,j) =  grid_.u(i+1 , j, WEST ) + ( grid_.dx()/ grid_.dy()) * ( grid_.v()( i , j ) - grid_.v( i , j-1 , NORTH) ) ;

		if (grid_.isEmpty(i-1,j-1) )
			grid_.v()(i-1,j-1) = grid_.v(i , j-1 , NORTH ) + ( grid_.dx()/ grid_.dy()) * ( grid_.u()( i , j ) - grid_.u( i , j-1 , NORTH) ) ;
		else ///////////////////////////////////////////// PROBLEM ////////////////////////////////////////////////////////////////////////
			grid_.v()(i-1,j-1) = grid_.v(i-1 , j-2 , NORTH ) - ( grid_.dy()/ grid_.dx()) * ( grid_.u( i-1 , j-1 , NORTH ) - grid_.u( i-2 , j-1 , NORTH) ) ;
	 }

	// NORTH : only the cell (i,j+1) is empty (North of cell(i,j) )

	else if ( grid_.isEmpty(i,j+1) )
	 {
		grid_.p()(i,j) =  ( 2 / (Re_*grid_.dy()) ) * ( grid_.v()( i , j) - grid_.v( i , j-1 , NORTH) ) ;
		grid_.v()(i,j) =  grid_.v(i , j-1 , NORTH ) - ( grid_.dy()/ grid_.dx()) * ( grid_.u()( i , j ) - grid_.u( i-1 , j , EAST) ) ;

		if (grid_.isEmpty(i-1,j+1) )
			grid_.u()(i-1,j+1) =  grid_.u(i-1 , j , EAST ) - ( grid_.dy()/ grid_.dx()) * ( grid_.v()( i , j ) - grid_.v( i-1 , j , EAST ) ) ;
		else ///////////////////////////////////////////// PROBLEM ////////////////////////////////////////////////////////////////////////
			grid_.u()(i-1,j+1) =  grid_.u(i-2 , j+1 , EAST ) - ( grid_.dx()/ grid_.dy()) * ( grid_.v( i-1 , j+1 , EAST ) - grid_.v( i-1 , j , EAST ) ) ;
	 }

	// SOUTH : only the cell (i,j-1) is empty (South of cell(i,j) )

	else if ( grid_.isEmpty(i,j-1) )
	 {
		grid_.p()(i , j) =  ( 2 / (Re_*grid_.dy()) ) * ( grid_.v()( i , j) - grid_.v( i , j-1 , NORTH) ) ;
		grid_.v()(i , j) =  grid_.v(i , j+1 , SOUTH ) + ( grid_.dy()/ grid_.dx()) * ( grid_.u()( i , j ) - grid_.u( i-1 , j , EAST) ) ;

		if (grid_.isEmpty(i-1,j-1) )
			grid_.u()(i-1,j-1) =  grid_.u(i-1 , j , EAST ) + ( grid_.dy()/ grid_.dx()) * ( grid_.v()( i , j ) - grid_.v( i-1 , j , EAST ) ) ;
		else ///////////////////////////////////////////// PROBLEM ////////////////////////////////////////////////////////////////////////
			grid_.u()(i-1,j-1) =  grid_.u(i-2 , j+1 , EAST ) - ( grid_.dx()/ grid_.dy()) * ( grid_.v( i-1 , j-1 , EAST ) - grid_.v( i-1 , j-2 , EAST ) ) ;
	 }
}
//*******************************************************************************************************************

void FluidSimulator::two_empty_neighbour(int i , int j , const int& dt)
{

}
//*******************************************************************************************************************

void FluidSimulator::three_empty_neighbour(int i , int j , const int& dt)
{

}
//*******************************************************************************************************************

void FluidSimulator::four_empty_neighbour(int i , int j , const int& dt)
{

}