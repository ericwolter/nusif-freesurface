
#include "FluidSimulator.hh"
#include <iostream>



// constructor for FluidSimulator
FluidSimulator::FluidSimulator(const FileReader &conf)
    : grid_(conf), solver_(conf)
{
    particle_tracer_ = ParticleTracer(&grid_);

    name_ = conf.getStringParameter("name");

    PROG("construct FluidSimulator with a file");
    gamma_ = conf.getRealParameter("gamma");
    CHECK_MSG((gamma_ >= 0 && gamma_ <= 1), "wrong input for gamma: " << gamma_);
    Re_ = conf.getRealParameter("Re");
    CHECK_MSG((Re_ > 0), "wrong input for Re: " << Re_);
    gx_ = conf.getRealParameter("GX");
    gy_ = conf.getRealParameter("GY");
    dt_ = conf.getRealParameter("dt");
    CHECK_MSG((dt_ > 0), "wrong input for dt: " << dt_);

    vel_N = conf.getRealParameter("boundary_velocity_N");
    vel_S = conf.getRealParameter("boundary_velocity_S");
    vel_E = conf.getRealParameter("boundary_velocity_E");
    vel_W = conf.getRealParameter("boundary_velocity_W");

    std::string no = conf.getStringParameter("boundary_condition_N");
    std::string so = conf.getStringParameter("boundary_condition_S");
    std::string we = conf.getStringParameter("boundary_condition_W");
    std::string es = conf.getStringParameter("boundary_condition_E");

    if (no == "inflow")
    {
        cond_N = INFLOW;
    }
    else if (no == "outflow")
    {
        ASSERT_MSG(vel_N == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_N = OUTFLOW;
    }
    else if (no == "free-slip")
    {
        cond_N = SLIP;
    }
    else
    {
        cond_N = NOSLIP;
    }

    if (so == "inflow")
    {
        cond_S = INFLOW;
    }
    else if (so == "outflow")
    {
        ASSERT_MSG(vel_S == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_S = OUTFLOW;
    }
    else if (so == "free-slip")
    {
        cond_S = SLIP;
    }
    else
    {
        cond_S = NOSLIP;
    }

    if (we == "inflow")
    {
        cond_W = INFLOW;
    }
    else if (we == "outflow")
    {
        ASSERT_MSG(vel_W == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_W = OUTFLOW;
    }
    else if (we == "free-slip")
    {
        cond_W = SLIP;
    }
    else
    {
        cond_W = NOSLIP;
    }

    if (es == "inflow")
    {
        cond_E = INFLOW;
    }
    else if (es == "outflow")
    {
        ASSERT_MSG(vel_E == 0.0, "there should be no velocity values for an outflow boundary!");
        cond_E = OUTFLOW;
    }
    else if (es == "free-slip")
    {
        cond_E = SLIP;
    }
    else
    {
        cond_E = NOSLIP;
    }

    safetyfac_ = conf.getRealParameter("safetyfactor");
    timeStepNr = (unsigned int)conf.getIntParameter("timesteps");
    CHECK_MSG((timeStepNr >= 0), "wrong input for timesteps: " << timeStepNr);
    outPutInt = (unsigned int)conf.getIntParameter("outputinterval");
    CHECK_MSG((outPutInt > 0), "wrong input for outputinterval: " << outPutInt);
    normfreq = (unsigned int)conf.getIntParameter("normalizationfrequency");
    CHECK_MSG((normfreq > 0), "wrong input for normalizationfrequency: " << normfreq);
    uInit_ = conf.getRealParameter("U_INIT");
    vInit_ = conf.getRealParameter("V_INIT");
    pInit_ = conf.getRealParameter("P_INIT");
    imax = conf.getIntParameter("imax");
    CHECK_MSG((imax > 0), "wrong input for imax: " << imax);
    jmax = conf.getIntParameter("jmax");
    CHECK_MSG((jmax > 0), "wrong input for jmax: " << jmax);

    // obstacles:
    // rectangle
    rectX_ = conf.getRealParameter("RectangleX1");
    CHECK_MSG((rectX_ >= 0), "wrong input for RectangleX1: " << rectX_);
    rectXX_ = conf.getRealParameter("RectangleX2");
    CHECK_MSG((rectXX_ >= 0), "wrong input for RectangleX2: " << rectXX_);
    rectY_ = conf.getRealParameter("RectangleY1");
    CHECK_MSG((rectY_ >= 0), "wrong input for RectangleY1: " << rectY_);
    rectYY_ = conf.getRealParameter("RectangleY2");
    CHECK_MSG((rectYY_ >= 0), "wrong input for RectangleY2: " << rectYY_);
    // circle
    circX_ = conf.getRealParameter("CircleX");
    CHECK_MSG((circX_ >= 0), "wrong input for CircleX: " << circX_);
    circY_ = conf.getRealParameter("CircleY");
    CHECK_MSG((circY_ >= 0), "wrong input for CircleY: " << circY_);
    circR_ = conf.getRealParameter("CircleR");
    CHECK_MSG((circR_ >= 0), "wrong input for CircleR: " << circR_);

    // particle
    // rectangle
    rectX1_particle_ = conf.getRealParameter("RectangleParticleX1");
    CHECK_MSG((rectX1_particle_ >= 0), "wrong input for RectangleParticleX1: " << rectX1_particle_);
    rectX2_particle_ = conf.getRealParameter("RectangleParticleX2");
    CHECK_MSG((rectX2_particle_ >= 0), "wrong input for RectangleParticleX2: " << rectX2_particle_);
    rectY1_particle_ = conf.getRealParameter("RectangleParticleY1");
    CHECK_MSG((rectY1_particle_ >= 0), "wrong input for RectangleParticleY1: " << rectY1_particle_);
    rectY2_particle_ = conf.getRealParameter("RectangleParticleY2");
    CHECK_MSG((rectY2_particle_ >= 0), "wrong input for RectangleParticleY2: " << rectY2_particle_);
    // circle
    circX_particle_ = conf.getRealParameter("CircleParticleX");
    CHECK_MSG((circX_particle_ >= 0), "wrong input for CircleParticleX: " << circX_particle_);
    circY_particle_ = conf.getRealParameter("CircleParticleY");
    CHECK_MSG((circY_particle_ >= 0), "wrong input for CircleParticleY: " << circY_particle_);
    circR_particle_ = conf.getRealParameter("CircleParticleR");
    CHECK_MSG((circR_particle_ >= 0), "wrong input for CircleParticleR: " << circR_particle_);

}

void FluidSimulator::composeRHS()
{
    // formula 3.38
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
            grid_.rhs()(i - 1, j - 1) = ((grid_.f()(i, j - 1) - grid_.f()(i - 1, j - 1)) / grid_.dx() + (grid_.g()(i - 1, j) - grid_.g()(i - 1, j - 1)) / grid_.dy()) / dt_;
    }
}

void FluidSimulator::updateVelocities()
{
    real dtdx = dt_ / grid_.dx();
    real dtdy = dt_ / grid_.dy();

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {

            if (grid_.isFluid(i, j))    // update obly for fluid cells
            {

                if (i != imax && grid_.isFluid(i + 1, j))   // formula 3.34
                    grid_.u()(i, j) = grid_.f()(i, j - 1) - (grid_.p(i + 1, j, WEST) - grid_.p()(i, j)) * dtdx;

                if (j != jmax && grid_.isFluid(i, j + 1))   // formula 3.35
                    grid_.v()(i, j) = grid_.g()(i - 1, j) - (grid_.p(i, j + 1, SOUTH) - grid_.p()(i, j)) * dtdy;

            }
        }
    }
}

void FluidSimulator::determineNextDT(real const &limit)
{
    // first stability condition
    real one = Re_ / (2 / pow(grid_.dx(), 2) + 2 / pow(grid_.dy(), 2));
    bool first = dt_ < one;
    // CFL / 3.49
    real two = grid_.dx() / grid_.u().maxE();
    bool second = dt_ < two;
    real three = grid_.dy() / grid_.v().maxE();
    bool third = dt_ < three;
    // control step size
    if (!(first && second && third))     // 3.50
    {
        real minimum = std::min(one, two);
        minimum = std::min(minimum, three);

        if (limit > 0 && limit <= 1)     // check safety factor
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
    if (cond_N == NOSLIP)     // no-slip
    {
        for (int i = 1; i <= imax; i++)
        {
            grid_.v()(i, jmax) = 0;
            grid_.u()(i, jmax + 1) = 2 * vel_N - grid_.u()(i, jmax);
        }

    }
    else if (cond_N == INFLOW)       // inflow
    {
        for (int i = 1; i <= imax; i++)
        {
            if (grid_.isFluid(i, jmax))
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
    else if (cond_N == SLIP)       // free-slip
    {
        for (int i = 1; i <= imax; i++)
        {
            if (grid_.isFluid(i, jmax))
            {
                grid_.v()(i, jmax) = vel_N;
                grid_.u()(i, jmax + 1) = grid_.u()(i, jmax);
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
        for (int i = 1; i <= imax; i++)
        {
            grid_.v()(i, jmax) = grid_.v()(i, jmax - 1);
            grid_.u()(i, jmax + 1) = grid_.u()(i, jmax);
        }
    }

    // conditions for south border
    if (cond_S == NOSLIP)     // no-slip
    {
        for (int i = 1; i <= imax; i++)
        {
            grid_.v()(i, 0) = 0;
            grid_.u()(i, 0) = 2 * vel_S - grid_.u()(i, 1);
        }

    }
    else if (cond_S == INFLOW)       // inflow
    {
        for (int i = 1; i <= imax; i++)
        {
            if (grid_.isFluid(i, 1))
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
    else if (cond_S == SLIP)        // free-slip
    {
        for (int i = 1; i <= imax; i++)
        {
            if (grid_.isFluid(i, 1))
            {
                grid_.v()(i, 0) = vel_S;
                grid_.u()(i, 0) = grid_.u()(i, 1);
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
        for (int i = 1; i <= imax; i++)
        {
            grid_.v()(i, 0) = grid_.v()(i, 1);
            grid_.u()(i, 0) = grid_.u()(i, 1);
        }
    }

    // conditions for west border
    if (cond_W == NOSLIP)     // no-slip
    {
        for (int j = 1; j <= jmax; j++)
        {
            grid_.u()(0, j) = 0;
            grid_.v()(0, j) = 2 * vel_W - grid_.v()(1, j);
        }

    }
    else if (cond_W == INFLOW)       // inflow
    {
        for (int j = 1; j <= jmax; j++)
        {
            if (grid_.isFluid(1, j))
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
    else if (cond_W == SLIP)       // free-slip
    {
        for (int j = 1; j <= jmax; j++)
        {
            if (grid_.isFluid(1, j))
            {
                grid_.u()(0, j) = vel_W;
                grid_.v()(0, j) = grid_.v()(1, j);
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
        for (int j = 1; j <= jmax; j++)
        {
            grid_.u()(0, j) = grid_.u()(1, j);
            grid_.v()(0, j) = grid_.v()(1, j);
        }
    }

    // conditions for east border
    if (cond_E == NOSLIP)     // no-slip
    {
        for (int j = 1; j <= jmax; j++)
        {
            grid_.u()(imax, j) = 0;
            grid_.v()(imax + 1, j) = 2 * vel_E - grid_.v()(imax, j);
        }

    }
    else if (cond_E == INFLOW)       // inflow
    {
        for (int j = 1; j <= jmax; j++)
        {
            if (grid_.isFluid(imax, j))
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
    else if (cond_E == SLIP)       // free-slip
    {
        for (int j = 1; j <= jmax; j++)
        {
            if (grid_.isFluid(imax, j))
            {
                grid_.u()(imax, j) = vel_E;
                grid_.v()(imax + 1, j) = grid_.v()(imax, j);
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
        for (int j = 1; j <= jmax; j++)
        {
            grid_.u()(imax, j) = grid_.u()(imax - 1, j);
            grid_.v()(imax + 1, j) = grid_.v()(imax, j);
        }
    }

}

void FluidSimulator::simulate(real duration)
{
    VTKWriter vtkWriter(name_);

    real t = 0;
    unsigned int n = 0;

    PROG("initialize u, v, p, rhs");
    if (name_ == "backstep")
    {
        int half = (int)(rectYY_ / grid_.dy());
        for (int i = 0; i < grid_.u().getSize(0); ++i)
        {
            for (int j = 0; j < grid_.u().getSize(1); ++j)
            {
                if (j > half)
                {
                    grid_.u()(i, j) = 1;
                }
                else
                {
                    grid_.u()(i, j) = 0;
                }
            }
        }
    }
    else
    {
        grid_.u().fill(uInit_);
    }
    grid_.v().fill(vInit_);
    grid_.p().fill(pInit_);
    grid_.rhs().fill(0);
    PROG("set inner obstacles");
    grid_.createRectangle(rectX_, rectY_, rectXX_, rectYY_);
    grid_.createCircle(circX_, circY_, circR_);
    PROG("set initial partciles");
    if (rectX1_particle_ + rectX2_particle_ + rectY1_particle_ + rectY2_particle_ + circR_particle_ + circX_particle_ + circY_particle_ == 0)
    {
        // fill non obstacle cells with particles
        for (int i = 1; i <= imax; ++i)
        {
            for (int j = 1; j <= jmax; ++j)
            {
                if (!grid_.isObstacle(i, j))
                    particle_tracer_.fillCell(i, j, grid_.ppc(), 0);
            }
        }
    }
    else
    {
        particle_tracer_.addRectangle((int)(rectX1_particle_ / grid_.dx()), (int)(rectX2_particle_ / grid_.dy()), (int)(rectY1_particle_ / grid_.dx()), (int)(rectY2_particle_ / grid_.dy()), 0);
        particle_tracer_.addCircle((int)(circX_particle_ / grid_.dx()), (int)(circY_particle_ / grid_.dy()), (int)(circR_particle_), 0);
    }

    PROG("initial refreshBoundaries");
    refreshBoundaries();

    while (t <= duration)
    {
        // without 0th step:
        // if ( n%outPutInt == 0 && n != 0 )
        if (n % outPutInt == 0)
            vtkWriter.write(grid_, &particle_tracer_);
        PROG(n << "'th timestep: determine next dt");
        determineNextDT(safetyfac_);
        //particle_tracer_.markCells();
        PROG(n << "'th timestep: set u, v, p at the free boundary");
        set_UVP_surface(dt_, true);
        computeFG();
        PROG(n << "'th timestep: compute the right-hand side of the pressure equation");
        composeRHS();
        PROG(n << "'th timestep: solve pressure equation");
        solv().solve(grid_);
        PROG(n << "'th timestep: update u and v in fluid domain");
        updateVelocities();
        PROG(n << "'th timestep: refresh boundaries");
        refreshBoundaries();
        PROG(n << "'th timestep: set u, v at the free boundary");
        set_UVP_surface(dt_, false);
        //particle_tracer_.advanceParticles(dt_);
        if (n % normfreq == 0)
            normalization();
        t += dt_;
        n++;
    }
}

void FluidSimulator::simulateTimeStepCount(unsigned int nrOfTimeSteps)
{
    VTKWriter vtkWriter(name_);
    unsigned int n = 0;

    PROG("set inner obstacles");
    //grid_.createRectangle(rectX_, rectY_, rectXX_, rectYY_);
    //grid_.createCircle(circX_, circY_, circR_);
    PROG("set initial particles");
    if (rectX1_particle_ + rectX2_particle_ + rectY1_particle_ + rectY2_particle_ + circR_particle_ + circX_particle_ + circY_particle_ == 0.0)
    {
        PROG("no explicit particles defined -> fill all");
        // fill non obstacle cells with particles
        for (int i = 1; i <= imax; ++i)
        {
            for (int j = 1; j <= jmax; ++j)
            {
                if (!grid_.isObstacle(i, j))
                {
                    particle_tracer_.fillCell(i, j, grid_.ppc(), 0);
                }
            }
        }
    }
    else
    {
        PROG("adding particles");
        particle_tracer_.addRectangle((int)(rectX1_particle_ / grid_.dx()), (int)(rectY1_particle_ / grid_.dy()), (int)(rectX2_particle_ / grid_.dx()), (int)(rectY2_particle_ / grid_.dy()), 0);
        particle_tracer_.addCircle((int)(circX_particle_ / grid_.dx()), (int)(circY_particle_ / grid_.dy()), (int)(circR_particle_), 0);
    }
    grid_.createPng("test.png");
    particle_tracer_.markCells();

    PROG("initialize u, v, p, rhs");
    if (name_ == "backstep")
    {
        int half = (int)(rectYY_ / grid_.dy());
        for (int i = 0; i < grid_.u().getSize(0); ++i)
        {
            for (int j = 0; j < grid_.u().getSize(1); ++j)
            {
                if (j > half && grid_.isFluid(i, j))
                {
                    grid_.u()(i, j) = 1.0;
                }
                else
                {
                    grid_.u()(i, j) = 0.0;
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < grid_.u().getSize(0); ++i)
        {
            for (int j = 0; j < grid_.u().getSize(1); ++j)
            {
                if (grid_.isFluid(i, j))
                {
                    grid_.u()(i, j) = uInit_ ;
                }
                else
                {
                    grid_.u()(i, j) = 0.0 ;
                }
            }
        }
    }

    for (int i = 0; i < grid_.v().getSize(0); ++i)
    {
        for (int j = 0; j < grid_.v().getSize(1); ++j)
        {
            if (grid_.isFluid(i, j))
            {
                grid_.v()(i, j) = vInit_ ;
            }
            else
            {
                grid_.v()(i, j) = 0.0 ;
            }
        }
    }

    for (int i = 0; i < grid_.p().getSize(0); ++i)
    {
        for (int j = 0; j < grid_.p().getSize(1); ++j)
        {
            if (grid_.isFluid(i, j))
            {
                grid_.p()(i, j) = pInit_ ;
            }
            else
            {
                grid_.p()(i, j) = 0.0 ;
            }
        }
    }

    grid_.rhs().fill(0);

    PROG("initial refreshBoundaries");
    refreshBoundaries();

    while (n <= nrOfTimeSteps)
    {
        // without 0th step:
        // if ( n%outPutInt == 0 && n != 0 )
        if (n % outPutInt == 0)
            vtkWriter.write(grid_, &particle_tracer_);
        PROG(n << "'th timestep: determine next dt");
        determineNextDT(safetyfac_);
        std::cout << dt_ << std::endl;
        particle_tracer_.markCells();
        grid_.obs().print();
        grid_.u().print();
        grid_.v().print();
        grid_.p().print();
        PROG(n << "'th timestep: set u, v, p at the free boundary");
        set_UVP_surface(dt_, true);
        grid_.u().print();
        grid_.v().print();
        grid_.p().print();
        computeFG();
        grid_.f().print();
        grid_.g().print();
        PROG(n << "'th timestep: compute the right-hand side of the pressure equation");
        composeRHS();
        PROG(n << "'th timestep: solve pressure equation");
        solv().solve(grid_);
        grid_.u().print();
        grid_.v().print();
        grid_.p().print();
        PROG(n << "'th timestep: update u and v in fluid domain");
        updateVelocities();
        grid_.u().print();
        grid_.v().print();
        grid_.p().print();
        PROG(n << "'th timestep: refresh boundaries");
        refreshBoundaries();
        grid_.u().print();
        grid_.v().print();
        grid_.p().print();
        PROG(n << "'th timestep: set u, v at the free boundary");
        set_UVP_surface(dt_, false);
        grid_.u().print();
        grid_.v().print();
        grid_.p().print();
        particle_tracer_.advanceParticles(dt_);
        if (n % normfreq == 0)
            normalization();
        n++;
    }

}

void FluidSimulator::normalization()
{
    real psum = 0;
    int count = 0;

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {

            if (grid_.isFluid(i, j))
            {
                psum += grid_.p()(i, j);
                count++;
            }
        }
    }

    real pMean = psum / grid_.getNumFluid();
    std::cout << "norm1: " << grid_.getNumFluid() << std::endl;
    std::cout << "norm2: " << count << std::endl;
    // substract the mean
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {

            if (grid_.isFluid(i, j))
                grid_.p()(i, j) -= pMean;
        }
    }

}

// help functions for derivatives
// d(u^2)/dx
real FluidSimulator::dxuu(int i, int j)
{
    ASSERT_MSG((i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG((i + 1 < grid_.u().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG((j >= 0), "wrong input for j: " << j);
    real res = (grid_.u()(i, j) + grid_.u(i + 1, j, WEST)) * (grid_.u()(i, j) + grid_.u(i + 1, j, WEST)) / (4 * grid_.dx());
    res -= (grid_.u(i - 1, j, EAST) + grid_.u()(i, j)) * (grid_.u(i - 1, j, EAST) + grid_.u()(i, j)) / (4 * grid_.dx());
    res += gamma_ * (fabs(grid_.u()(i, j) + grid_.u(i + 1, j, WEST)) * (grid_.u()(i, j) - grid_.u(i + 1, j, WEST)) * 0.25) / grid_.dx();
    res -= gamma_ * (fabs(grid_.u(i - 1, j, EAST) + grid_.u()(i, j)) * (grid_.u(i - 1, j, EAST) - grid_.u()(i, j)) * 0.25) / grid_.dx();

    return res;
}

// d(v^2)/dy
real FluidSimulator::dyvv(int i, int j)
{
    ASSERT_MSG((i >= 0), "wrong input for i: " << i);
    ASSERT_MSG((j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG((j + 1 < grid_.v().getSize(1)), "wrong input for j: " << j);
    real res = (grid_.v()(i, j) + grid_.v(i, j + 1, SOUTH)) * (grid_.v()(i, j) + grid_.v(i, j + 1, SOUTH)) / (4 * grid_.dy());
    res -= (grid_.v(i, j - 1, NORTH) + grid_.v()(i, j)) * (grid_.v(i, j - 1, NORTH) + grid_.v()(i, j)) / (4 * grid_.dy());
    res += gamma_ * (fabs(grid_.v()(i, j) + grid_.v(i, j + 1, SOUTH)) * (grid_.v()(i, j) - grid_.v(i, j + 1, SOUTH)) * 0.25) / grid_.dy();
    res -= gamma_ * (fabs(grid_.v(i, j - 1, NORTH) + grid_.v()(i, j)) * (grid_.v(i, j - 1, NORTH) - grid_.v()(i, j)) * 0.25) / grid_.dy();

    return res;
}
////////////////////////////////////////////////////////////////////////////////
// d(uv)/dy
real FluidSimulator::dyuv(int i, int j)
{
    ASSERT_MSG((i >= 0), "wrong input for i: " << i);
    ASSERT_MSG((i + 1 < grid_.v().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG((j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG((j + 1 < grid_.u().getSize(1)), "wrong input for j: " << j);
    real diag = grid_.v(i + 1, j - 1, DIAG);
    //     real diag = 0.0;
    //     if (grid_.isFluid(i + 1, j - 1))
    //     {
    //         diag = grid_.v()(i + 1, j - 1);
    //     }
    //     else
    //     {
    //         diag = 0.5 * (grid_.v(i + 1, j, WEST) + grid_.v(i, j - 1, NORTH));
    //     }
    real res = ((grid_.v()(i, j) + grid_.v(i + 1, j, WEST)) * (grid_.u()(i, j) + grid_.u(i, j + 1, SOUTH)) * 0.25) / grid_.dy();
    res -= ((grid_.v(i, j - 1, NORTH) + diag) * (grid_.u(i, j - 1, NORTH) + grid_.u()(i, j)) * 0.25) / grid_.dy();
    res += gamma_ * (fabs(grid_.v()(i, j) + grid_.v(i + 1, j, WEST)) * (grid_.u()(i, j) - grid_.u(i, j + 1, SOUTH)) * 0.25) / grid_.dy();
    res -= gamma_ * (fabs(grid_.v(i, j - 1, NORTH) + diag) * (grid_.u(i, j - 1, NORTH) - grid_.u()(i, j)) * 0.25) / grid_.dy();

    return res;
}

// d(uv)/dx
real FluidSimulator::dxuv(int i, int j)
{
    ASSERT_MSG((j >= 0), "wrong input for j: " << j);
    ASSERT_MSG((j + 1 < grid_.u().getSize(1)), "wrong input for j: " << j);
    ASSERT_MSG((i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG((i + 1 < grid_.v().getSize(0)), "wrong input for i: " << i);
    real diag = grid_.u(i - 1, j + 1, DIAG);
    //     real diag = 0.0;
    //     if (grid_.isFluid(i - 1, j + 1))
    //     {
    //         diag = grid_.u()(i - 1, j + 1);
    //     }
    //     else
    //     {
    //         diag = 0.5 * (grid_.u(i, j + 1, EAST) + grid_.u(i - 1, j, SOUTH));
    //     }
    real res = ((grid_.u()(i, j) + grid_.u(i, j + 1, SOUTH)) * (grid_.v()(i, j) + grid_.v(i + 1, j, WEST)) * 0.25) / grid_.dx();
    res -= ((grid_.u(i - 1, j, EAST) + diag) * (grid_.v(i - 1, j, EAST) + grid_.v()(i, j)) * 0.25) / grid_.dx();
    res += gamma_ * (fabs(grid_.u()(i, j) + grid_.u(i, j + 1, SOUTH)) * (grid_.v()(i, j) - grid_.v(i + 1, j, WEST)) * 0.25) / grid_.dx();
    res -= gamma_ * (fabs(grid_.u(i - 1, j, EAST) + diag) * (grid_.v(i - 1, j, EAST) - grid_.v()(i, j)) * 0.25) / grid_.dx();

    return res;
}

// d^2u/dx^2
real FluidSimulator::ddxu(int i, int j)
{
    ASSERT_MSG((i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG((i + 1 < grid_.u().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG((j >= 0), "wrong input for j: " << j);
    real res = (grid_.u(i + 1, j, WEST) - 2 * grid_.u()(i, j) + grid_.u(i - 1, j, EAST)) / (grid_.dx() * grid_.dx());

    return res;
}

// d^2v/dx^2
real FluidSimulator::ddxv(int i, int j)
{
    ASSERT_MSG((i - 1 >= 0), "wrong input for i: " << i);
    ASSERT_MSG((i + 1 < grid_.v().getSize(0)), "wrong input for i: " << i);
    ASSERT_MSG((j >= 0), "wrong input for j: " << j);
    real res = (grid_.v(i + 1, j, WEST) - 2 * grid_.v()(i, j) + grid_.v(i - 1, j, EAST)) / (grid_.dx() * grid_.dx());

    return res;
}

// d^2u/dy^2
real FluidSimulator::ddyu(int i, int j)
{
    ASSERT_MSG((j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG((j + 1 < grid_.u().getSize(1)), "wrong input for j: " << j);
    ASSERT_MSG((i >= 0), "wrong input for i: " << i);
    real res = (grid_.u(i, j + 1, SOUTH) - 2 * grid_.u()(i, j) + grid_.u(i, j - 1, NORTH)) / (grid_.dy() * grid_.dy());

    return res;
}

// d^2v/dy^2
real FluidSimulator::ddyv(int i, int j)
{
    ASSERT_MSG((j - 1 >= 0), "wrong input for j: " << j);
    ASSERT_MSG((j + 1 < grid_.v().getSize(1)), "wrong input for j: " << j);
    ASSERT_MSG((i >= 0), "wrong input for i: " << i);
    real res = (grid_.v(i, j + 1, SOUTH) - 2 * grid_.v()(i, j) + grid_.v(i, j - 1, NORTH)) / (grid_.dy() * grid_.dy());

    return res;
}

void FluidSimulator::computeFG()
{
    grid_.f().fill(0);
    grid_.g().fill(0);
    PROG("set boundary values for F");
    // boundary values for F
    for (int j = 1; j <= jmax; j++)
    {
        grid_.f()(0, j - 1) = grid_.u()(0, j);
        grid_.f()(imax, j - 1) = grid_.u()(imax, j);
    }

    PROG("set boundary values for G");
    // boundary values for G
    for (int i = 1; i <= imax; i++)
    {
        grid_.g()(i - 1, 0) = grid_.v()(i, 0);
        grid_.g()(i - 1, jmax) = grid_.v()(i, jmax);
    }

    real dtRe = dt_ / Re_;
    real dtGx = dt_ * gx_;
    real dtGy = dt_ * gy_;

    PROG("compute F");
    // compute F
    for (int i = 1; i < imax; i++)
    {
        for (int j = 1; j <= jmax; j++)
        {

            if (grid_.isFluid(i, j) && grid_.isFluid(i + 1, j))
                grid_.f()(i, j - 1) = grid_.u()(i, j) + dtRe * (ddxu(i, j) + ddyu(i, j)) - dt_ * (dxuu(i, j) + dyuv(i, j)) + dtGx;

        }
    }

    PROG("compute G");
    // compute G
    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j < jmax; j++)
        {

            if (grid_.isFluid(i, j) && grid_.isFluid(i, j + 1))
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

void FluidSimulator::set_UVP_surface(const real &dt, bool compP)
{

    for (int i = 1; i <= imax; i++)
    {
        for (int j = 1; j < jmax; j++)   // !?!
        {

            set_UVP_surface(i, j, dt, compP);
        }
    }
}


void FluidSimulator::set_UVP_surface(int i, int j , const real &dt, bool compP)
{

    int number_of_empty_neighbour = 0 ;

    if (grid_.isFluid(i, j))
    {
        // Count the number of empty neighbour cell
        if (grid_.isEmpty(i + 1, j))
            ++ number_of_empty_neighbour ;
        if (grid_.isEmpty(i - 1, j))
            ++ number_of_empty_neighbour ;
        if (grid_.isEmpty(i, j + 1))
            ++ number_of_empty_neighbour ;
        if (grid_.isEmpty(i, j - 1))
            ++ number_of_empty_neighbour ;

    }
    switch (number_of_empty_neighbour)
    {
    case (1):
        one_empty_neighbour(i , j , dt, compP) ;
        break;
    case (2):
        two_empty_neighbour(i , j , dt, compP) ;
        break;
    case (3):
        three_empty_neighbour(i , j , dt, compP) ;
        break;
    case (4):
        four_empty_neighbour(i , j , dt, compP) ;
        break;
    }
}

// real FluidSimulator::set_U_surface(int i, int j , const real &dt) // not finished yet!!!
// {
//     int number_of_empty_neighbour = 0 ;
//
//     if ( grid_.isFluid(i, j) )
//     {
//         // Count the number of empty neighbour cell
//         if ( grid_.isEmpty(i + 1, j) )
//             ++ number_of_empty_neighbour ;
//         if ( grid_.isEmpty(i - 1, j) )
//             ++ number_of_empty_neighbour ;
//         if ( grid_.isEmpty(i, j + 1) )
//             ++ number_of_empty_neighbour ;
//         if ( grid_.isEmpty(i, j - 1) )
//             ++ number_of_empty_neighbour ;
//
//     }
//
//     switch (number_of_empty_neighbour )
//     {
//
//     case (1): // one empty neighbour in the
//       if ( grid_.isEmpty(i+1, j) ) // east
//  return  grid_.u()(i - 1 , j) - ( grid_.dx() / grid_.dy()) * ( grid_.v()( i , j ) - grid_.v()( i , j - 1 ) ) ;
//       else // south/north/west
//  return grid_.u()(i,j) ;
//
//     case (2): // two empty neighbours in the
//       if (grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j - 1) ) // east and south
//         return grid_.u()( i - 1 , j ) ;
//       else if (grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j + 1) ) // east and north
//         return grid_.u()(i - 1 , j) ;
//       else if (grid_.isEmpty(i + 1, j) && grid_.isEmpty(i - 1, j)) // east and west
//         return grid_.u()(i, j) + gx_ * dt ;
//       else // north and south
//  return grid_.u()(i,j) ;
//
//     case (3): // three empty neighbours in the
//       if ( grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j + 1) && grid_.isEmpty(i - 1, j)) // east, north and west
//         return grid_.u()(i, j) + gx_ * dt ;
//       else if ( grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j - 1) && grid_.isEmpty(i - 1, j)) // east, south and west
//         return grid_.u()(i, j) + gx_ * dt ;
//       else if ( grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j - 1) && grid_.isEmpty(i, j + 1) ) // east, south and north
//         return grid_.u()(i - 1, j) - (grid_.dx() / grid_.dy() ) * ( grid_.v()(i, j) - grid_.v()(i, j - 1) ) ;
//       else // north, south and west
//  return grid_.u()(i,j) ;
//
//     case (4): // four empty neighbours
//       return grid_.u()(i, j) + gx_ * dt ;
//
//     }
// }
//
// real FluidSimulator::set_V_surface(int i, int j , const real &dt) already finished!
// {
//     int number_of_empty_neighbour = 0 ;
//
//     if ( grid_.isFluid(i, j) )
//     {
//         // Count the number of empty neighbour cell
//         if ( grid_.isEmpty(i + 1, j) )
//             ++ number_of_empty_neighbour ;
//         if ( grid_.isEmpty(i - 1, j) )
//             ++ number_of_empty_neighbour ;
//         if ( grid_.isEmpty(i, j + 1) )
//             ++ number_of_empty_neighbour ;
//         if ( grid_.isEmpty(i, j - 1) )
//             ++ number_of_empty_neighbour ;
//
//     }
//
//     switch (number_of_empty_neighbour )
//     {
//
//     case (1): // one empty neighbour in the
//       if ( grid_.isEmpty(i, j + 1) ) // north
//  return  grid_.v()(i , j - 1) - ( grid_.dy() / grid_.dx()) * ( grid_.u()( i , j ) - grid_.u()( i - 1 , j ) ) ;
//       else // south/east/west
//  return grid_.v()(i,j) ;
//
//     case (2): // two empty neighbours in the
//       if (grid_.isEmpty(i - 1, j) && grid_.isEmpty(i, j + 1) ) // west and north
//  return grid_.v()( i  , j - 1 ) ;
//       else if (grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j + 1) ) // east and north
//         return grid_.v()( i  , j - 1 ) ;
//       else if (grid_.isEmpty(i, j + 1) && grid_.isEmpty(i, j - 1)) // south and north
//         return grid_.v()(i, j) + gy_ * dt ;
//       else // east and west
//  return grid_.v()(i,j) ;
//
//     case (3): // three empty neighbours in the
//       if ( grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j + 1)  && grid_.isEmpty(i - 1, j)) // east, north and west
//         return grid_.v()(i, j) = grid_.v()(i, j - 1) - (grid_.dy() / grid_.dx() ) * ( grid_.u()(i, j) - grid_.u()(i - 1, j) ) ;
//       else if ( grid_.isEmpty(i - 1, j) && grid_.isEmpty(i, j - 1)  && grid_.isEmpty(i, j + 1) ) // west, south and east
//         return grid_.v()(i, j) + gy_ * dt ;
//       else if ( grid_.isEmpty(i + 1, j) && grid_.isEmpty(i, j - 1)  && grid_.isEmpty(i, j + 1) ) // east, south and north
//  return grid_.v()(i, j) + gy_ * dt ;
//       else // east, south and west
//  return grid_.v()(i,j) ;
//
//     case (4): // four empty neighbours
//       return grid_.v()(i, j) + gy_ * dt ;
//
//     }
// }


//*******************************************************************************************************************

void FluidSimulator::one_empty_neighbour(int i , int j , const real &dt, bool compP)
{
    // According to page 92 ,case 1 of the book
    // If cell (i,j) is Fluid, one of the other neighbour(North,East,West,South) is empty cell.
    // In any case of empty neighbour cell the formula is the same but we need to change for indices for mentioning
    // right neighbour cell.

    // EAST : only the cell (i+1,j) is empty (exactly like one which has been mentioned in the book)

    if (grid_.isEmpty(i + 1, j))
    {
        grid_.u()(i, j) =  grid_.u(i - 1 , j, EAST) - (grid_.dx() / grid_.dy()) * (grid_.v()(i , j) - grid_.v(i , j - 1, NORTH)) ;        
		if (compP)
            grid_.p()(i, j) = (2 / (Re_ * grid_.dx())) * (grid_.u()(i , j) - grid_.u(i - 1 , j, EAST)) ;

        if (grid_.isEmpty(i + 1, j - 1))
            grid_.v()(i + 1, j - 1) = grid_.v(i , j - 1, NORTH) - (grid_.dx() / grid_.dy()) * (grid_.u()(i , j) - grid_.u(i , j - 1, NORTH)) ;
        //         else ///////////////////////////////////////////// conti ////////////////////////////////////////////////////////////////////////
        //             grid_.v()(i + 1, j - 1) = grid_.v(i + 1 , j - 2, DIAG) - ( grid_.dy() / grid_.dx()) * ( grid_.u( i + 1 , j - 1, DIAG) - grid_.u( i , j - 1, NORTH ) ) ;

    }

    // WEST : only the cell (i-1,j) is empty (West of cell(i,j) )

    else if (grid_.isEmpty(i - 1, j))
    {
        grid_.u()(i - 1, j) =  grid_.u()(i , j) + (grid_.dx() / grid_.dy()) * (grid_.v()(i , j) - grid_.v(i , j - 1, NORTH)) ;        
		if (compP)
            grid_.p()(i, j)   = (2 / (Re_ * grid_.dx())) * (grid_.u()(i , j) - grid_.u(i - 1, j, EAST)) ;

        if (grid_.isEmpty(i - 1, j - 1))
            grid_.v()(i - 1, j - 1) = grid_.v(i, j - 1, NORTH) + (grid_.dx() / grid_.dy()) * (grid_.u()(i , j) - grid_.u(i, j - 1, NORTH)) ;
        //         else ///////////////////////////////////////////// conti ////////////////////////////////////////////////////////////////////////
        //             grid_.v()(i - 1, j - 1) = grid_.v(i - 1 , j - 2, DIAG) - ( grid_.dy() / grid_.dx()) * ( grid_.u( i - 1 , j - 1, DIAG ) - grid_.u( i - 2 , j - 1, DIAG ) ) ;
    }

    // NORTH : only the cell (i,j+1) is empty (North of cell(i,j) )

    else if (grid_.isEmpty(i, j + 1))
    {
		grid_.v()(i, j) =  grid_.v(i, j - 1, NORTH) - (grid_.dy() / grid_.dx()) * (grid_.u()(i , j) - grid_.u(i - 1, j, EAST)) ;
        if (compP)
            grid_.p()(i, j) = (2 / (Re_ * grid_.dy())) * (grid_.v()(i , j) - grid_.v(i, j - 1, NORTH)) ;

        if (grid_.isEmpty(i - 1, j + 1))
            grid_.u()(i - 1, j + 1) =  grid_.u(i - 1, j, EAST) - (grid_.dy() / grid_.dx()) * (grid_.v()(i , j) - grid_.v(i - 1, j, EAST)) ;
        //         else ///////////////////////////////////////////// conti + no values(?) ////////////////////////////////////////////////////////////
        //             grid_.u()(i - 1, j + 1) =  grid_.u(i - 2 , j + 1, DIAG) - ( grid_.dx() / grid_.dy()) * ( grid_.v( i - 1 , j + 1, DIAG ) - grid_.v( i - 1 , j, EAST ) ) ;
    }

    // SOUTH : only the cell (i,j-1) is empty (South of cell(i,j) )

    else if (grid_.isEmpty(i, j - 1))
    {
		grid_.v()(i, j - 1) =  grid_.v()(i , j) + (grid_.dy() / grid_.dx()) * (grid_.u()(i , j) - grid_.u(i - 1, j, EAST)) ;
        if (compP)
            grid_.p()(i, j) = (2 / (Re_ * grid_.dy())) * (grid_.v()(i , j) - grid_.v(i, j - 1, NORTH)) ;    

        if (grid_.isEmpty(i - 1, j - 1))
            grid_.u()(i - 1, j - 1) =  grid_.u(i - 1, j, EAST) + (grid_.dy() / grid_.dx()) * (grid_.v()(i , j) - grid_.v(i - 1, j, EAST)) ;
        //         else ///////////////////////////////////////////// conti + no values(?) //////////////////////////////////////////////////////////////
        //             grid_.u()(i - 1, j - 1) =  grid_.u(i - 2 , j + 1, DIAG ) - ( grid_.dx() / grid_.dy()) * ( grid_.v( i - 1 , j - 1, DIAG ) - grid_.v( i - 1 , j - 2, DIAG ) ) ;
    }
}
//*******************************************************************************************************************

void FluidSimulator::two_empty_neighbour(int i , int j , const real &dt, bool compP)
{

    real Re_inverse = (real) 1.0 / Re_ ;
    real dy_inverse = (real) 1.0 / grid_.dy() ;
    real dx_inverse = (real) 1.0 / grid_.dx() ;

    // According to page 93 and 95 ,this case divides into two parts :
    // 1- two empty cells have a corner in share , e.g. Cell(i+1,j) and Cell(i,j-1) are two empty cell
    // 2- two empty cells are on the opposite, e.g. Cell(i+1,j) and Cell(i-1,j)

    // Case 1 : in this case we have four options for our two empty neighbour

    // Cell(i+1,j) and Cell(i,j-1) are empty cell (Like one in book page 94 ) SE

    if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i, j - 1))
    {
        grid_.u()(i, j)   = grid_.u(i - 1, j, EAST) ;
        grid_.v()(i, j - 1)   = grid_.v()(i   , j);
        if (compP)
            grid_.p()(i, j)   = -0.5 * Re_inverse * (
                                    (grid_.u(i, j + 1, SOUTH) + grid_.u(i - 1, j + 1, DIAG) - grid_.u()(i, j) - grid_.u(i - 1, j, EAST)) * dy_inverse
                                    +
                                    (grid_.v()(i, j) + grid_.v(i, j - 1, NORTH) - grid_.v(i - 1, j, EAST) - grid_.v(i - 1, j - 1, DIAG)) * dx_inverse
                                );
        if (grid_.isEmpty(i - 1, j - 1))
            grid_.u()(i - 1, j - 1) = grid_.u(i - 1, j, EAST) + (grid_.dy() / grid_.dx()) * (grid_.v(i, j - 1, NORTH) - grid_.v(i - 1, j - 1, DIAG))  ;
        if (grid_.isEmpty(i + 1, j + 1))  /////////////////////////// no values ////////////////////////////////
            grid_.v()(i + 1, j)   = grid_.v()(i, j) - (grid_.dx() / grid_.dy()) * (grid_.u(i, j + 1, SOUTH) - grid_.u()(i, j)) ;
        grid_.u()(i, j - 1)   = grid_.u()(i, j)     ;
        grid_.v()(i + 1, j - 1) = grid_.v(i, j - 1, NORTH) ;
    }

    // Cell(i-1,j) and Cell(i,j+1) are empty cell  NW
    else if (grid_.isEmpty(i - 1, j)  && grid_.isEmpty(i, j + 1))
    {
        grid_.u()(i - 1, j)   = grid_.u()(i , j) ;
        grid_.v()(i, j)     = grid_.v()(i  , j - 1) ;
        if (compP)
            grid_.p()(i, j)   = -0.5 * Re_inverse * (
                                    (grid_.u()(i, j) + grid_.u(i - 1, j, EAST) - grid_.u(i, j - 1, NORTH) - grid_.u(i - 1, j - 1, DIAG)) * dy_inverse
                                    +
                                    (grid_.v(i + 1, j, WEST) - grid_.v(i + 1, j - 1, DIAG) - grid_.v()(i, j) - grid_.v(i, j - 1, NORTH)) * dx_inverse
                                ) ;
        if (grid_.isEmpty(i - 1, j - 1))
            grid_.v()(i - 1, j - 1) = grid_.v(i, j - 1, NORTH) + (grid_.dx() / grid_.dy()) * (grid_.u(i - 1, j, EAST) - grid_.u(i - 1 , j - 1, DIAG))  ;
        if (grid_.isEmpty(i + 1, j + 1))  ///////////////// no values ///////////////
            grid_.u()(i, j + 1)   = grid_.u()(i, j) - (grid_.dy() / grid_.dx()) * (grid_.v(i + 1, j, WEST) - grid_.v()(i, j)) ;
        grid_.u()(i - 1, j + 1)   = grid_.u(i - 1, j, EAST) ;
        grid_.v()(i - 1, j)     = grid_.v()(i, j) ;
    }

    // Cell(i,j+1) and Cell(i+1,j) are empty cell NE
    else if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i, j + 1))
    {
        grid_.u()(i, j)  = grid_.u(i - 1, j, EAST) ;
        grid_.v()(i, j)  = grid_.v(i, j - 1, NORTH) ;
        if (compP)
            grid_.p()(i, j)   =  0.5 * Re_inverse * (
                                     (grid_.u()(i, j) + grid_.u(i, j - 1, NORTH) - grid_.u(i - 1, j - 1, DIAG) - grid_.u(i - 1, j, EAST)) * dy_inverse
                                     +
                                     (grid_.v()(i, j) + grid_.v(i, j - 1, NORTH) - grid_.v(i - 1, j, EAST) - grid_.v(i - 1, j - 1, DIAG)) * dx_inverse
                                 ) ;
        if (grid_.isEmpty(i + 1, j - 1))
            grid_.v()(i + 1, j - 1) = grid_.v(i, j - 1, NORTH) - (grid_.dx() / grid_.dy()) * (grid_.u()(i , j) - grid_.u(i, j - 1, NORTH))  ;
        if (grid_.isEmpty(i - 1, j + 1))
            grid_.u()(i - 1, j + 1) = grid_.u(i - 1, j, EAST) - (grid_.dy() / grid_.dx()) * (grid_.v()(i, j) - grid_.v(i - 1, j, EAST)) ;
        grid_.u()(i, j + 1)   = grid_.u()(i, j) ;
        grid_.v()(i + 1, j)   = grid_.v()(i, j) ;
    }

    // Cell(i,j-1) and Cell(i-1,j) are empty cell SW
    else if (grid_.isEmpty(i - 1, j)  && grid_.isEmpty(i, j - 1))
    {
        grid_.u()(i - 1, j)  = grid_.u()(i , j) ;
        grid_.v()(i, j - 1)  = grid_.v()(i , j) ;
        if (compP)
            grid_.p()(i, j)   =  0.5 * Re_inverse * (
                                     (grid_.u(i - 1, j + 1, DIAG) + grid_.u(i, j + 1, SOUTH) - grid_.u(i - 1, j, EAST) - grid_.u()(i, j)) * dy_inverse
                                     +
                                     (grid_.v(i + 1, j, WEST) + grid_.v(i + 1, j - 1, DIAG) - grid_.v()(i, j) - grid_.v(i, j - 1, NORTH)) * dx_inverse
                                 ) ;
        if (grid_.isEmpty(i + 1, j - 1))
            grid_.u()(i, j - 1) = grid_.u()(i , j) + (grid_.dy() / grid_.dx()) * (grid_.v(i + 1, j - 1, DIAG) - grid_.v(i, j - 1, NORTH))  ;
        if (grid_.isEmpty(i - 1, j + 1))  ////////////////////// no values //////////////////////////
            grid_.v()(i - 1, j) = grid_.v()(i, j) + (grid_.dx() / grid_.dy()) * (grid_.u(i - 1, j + 1, DIAG) - grid_.u(i - 1, j, EAST)) ;
        grid_.u()(i - 1, j - 1)   = grid_.u(i - 1, j, EAST) ;
        grid_.v()(i - 1, j - 1)   = grid_.v(i, j - 1, NORTH) ;
    }

    // Case 2 : in this case we have two options for our two empty opposite neighbour

    // Cell(i+1,j) and Cell(i-1,j) are empty cell WE
    else if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i - 1, j))
    {
        grid_.u()(i, j) += gx_ * dt ;
        grid_.u()(i - 1, j) = grid_.u(i - 1, j, EAST) + gx_ * dt ;
        if (grid_.isEmpty(i - 1, j - 1))
            grid_.v()(i - 1 , j - 1) = grid_.v(i, j - 1, NORTH) + (grid_.dx() / grid_.dy()) * (grid_.u(i - 1, j, EAST) - grid_.u()(i - 1, j - 1)) ;
        if (grid_.isEmpty(i + 1, j - 1))
            grid_.v()(i + 1 , j - 1) = grid_.v(i, j - 1, NORTH) - (grid_.dx() / grid_.dy()) * (grid_.u()(i, j) - grid_.u(i, j - 1, NORTH)) ;
        if (compP)
            grid_.p()(i, j) = 0.0 ;
    }
    // Cell(i,j+1) and Cell(i,j-1) are empty cell NS
    else if (grid_.isEmpty(i, j + 1)  && grid_.isEmpty(i, j - 1))
    {
        grid_.v()(i, j) += gy_ * dt ;
        grid_.v()(i, j - 1) = grid_.v(i, j - 1, NORTH) + gy_ * dt ;
        if (grid_.isEmpty(i - 1, j - 1))
            grid_.u()(i - 1 , j - 1) = grid_.u(i - 1, j, EAST) + (grid_.dy() / grid_.dx()) * (grid_.v(i, j - 1, NORTH) - grid_.v(i - 1, j - 1, DIAG)) ;
        if (grid_.isEmpty(i - 1, j + 1))
            grid_.u()(i - 1 , j + 1) = grid_.u(i - 1, j, EAST) - (grid_.dy() / grid_.dx()) * (grid_.v()(i, j) - grid_.v(i - 1, j, EAST)) ;
        if (compP)
            grid_.p()(i, j) = 0.0 ;
    }

}
//*******************************************************************************************************************

void FluidSimulator::three_empty_neighbour(int i , int j , const real &dt, bool compP)
{
    std::cout << "three_empty_neighbour: " << i << ", " << j << std::endl;
    // According to page 96 of book
    // Generally we have four option for our three empty neighbours

    // cell(i+1,j) , cell(i,j+1) , cell(i-1,j) are empty WNE
    if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i, j + 1)  && grid_.isEmpty(i - 1, j))
    {
        std::cout << "WNE" << std::endl;
        if (compP)
            grid_.p()(i, j) = 0.0 ;
        grid_.u()(i, j) += gx_ * dt ;
        grid_.u()(i - 1, j) = grid_.u(i - 1, j, EAST) + gx_ * dt ;
        grid_.v()(i, j) = grid_.v(i, j - 1, NORTH) - (grid_.dy() / grid_.dx()) * (grid_.u()(i, j) - grid_.u(i - 1, j, EAST)) ;

		if (grid_.isEmpty(i - 1, j + 1))
		{
            grid_.u()(i - 1, j + 1) = grid_.u(i - 1, j, EAST) ;
			grid_.v()(i - 1, j)   = grid_.v()(i, j) ;
		}
		if (grid_.isEmpty(i + 1, j + 1))
		{
            grid_.u()(i, j + 1)   = grid_.u()(i, j) ;
			grid_.v()(i + 1, j)   = grid_.v()(i, j) ;
		}

        if (grid_.isEmpty(i + 1, j - 1))
            grid_.v()(i + 1, j - 1) = grid_.v(i, j - 1, NORTH) - (grid_.dx() / grid_.dy()) * (grid_.u()(i, j) - grid_.u(i, j - 1, NORTH)) ;

        if (grid_.isEmpty(i - 1, j - 1))
            grid_.v()(i - 1, j - 1) = grid_.v(i, j - 1, NORTH) + (grid_.dx() / grid_.dy()) * (grid_.u(i - 1, j, EAST) - grid_.u()(i - 1, j - 1)) ;
    }

    // cell(i+1,j) , cell(i,j-1) , cell(i-1,j) are empty WSE
    if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i, j - 1)  && grid_.isEmpty(i - 1, j))
    {
        std::cout << "WSE" << std::endl;
        std::cout << "grid_.v()(i, j): " << grid_.v()(i, j) << std::endl;
        std::cout << "grid_.u()(i, j): " << grid_.u()(i, j) << std::endl;
        std::cout << "grid_.u()(i - 1, j): " << grid_.u()(i - 1, j) << std::endl;
        std::cout << "grid_.v()(i - 1, j): " << grid_.v()(i - 1, j) << std::endl;

        if (compP)
            grid_.p()(i, j) = 0.0 ;
        grid_.u()(i, j) += gx_ * dt ;
        grid_.u()(i - 1, j) = grid_.u(i - 1, j, EAST) + gx_ * dt ;
        grid_.v()(i, j - 1) = grid_.v()(i, j) + (grid_.dy() / grid_.dx()) * (grid_.u()(i, j) - grid_.u(i - 1, j, EAST)) ;

		if (grid_.isEmpty(i - 1, j - 1))
		{
			grid_.u()(i - 1, j - 1)   = grid_.u(i - 1, j, EAST) ;
			grid_.v()(i - 1, j - 1)   = grid_.v(i, j - 1, NORTH) ;
		}
		if (grid_.isEmpty(i + 1, j - 1))
		{
			grid_.u()(i, j - 1)     = grid_.u()(i, j) ;
			grid_.v()(i + 1, j - 1)   = grid_.v(i, j - 1, NORTH) ;
		}
        

//         if (grid_.isEmpty(i + 1, j + 1))   /////////////// no values ////////////////////
//             grid_.v()(i + 1, j) = grid_.v()(i, j) + (grid_.dx() / grid_.dy()) * (grid_.u(i, j + 1, SOUTH) - grid_.u()(i, j)) ;
// 
//         if (grid_.isEmpty(i - 1, j + 1))
//             grid_.v()(i - 1, j) = grid_.v()(i, j) - (grid_.dx() / grid_.dy()) * (grid_.u(i - 1, j + 1, DIAG) - grid_.u(i - 1, j, EAST)) ;
    }

    // cell(i-1,j) , cell(i,j-1) , cell(i,j+1) are empty WSN
    if (grid_.isEmpty(i - 1, j)  && grid_.isEmpty(i, j - 1)  && grid_.isEmpty(i, j + 1))
    {
        if (compP)
            grid_.p()(i, j) = 0.0 ;
        grid_.v()(i, j) += gy_ * dt ;
        grid_.v()(i, j - 1) = grid_.v(i, j - 1, NORTH) + gy_ * dt ;
        grid_.u()(i - 1, j) = grid_.u()(i, j) + (grid_.dx() / grid_.dy()) * (grid_.v()(i, j) - grid_.v(i, j - 1, NORTH)) ;

		if (grid_.isEmpty(i - 1, j + 1))
		{
			grid_.u()(i - 1, j + 1)   = grid_.u(i - 1, j, EAST) ;
			grid_.v()(i - 1, j)     = grid_.v()(i, j) ;
		}
		if (grid_.isEmpty(i - 1, j - 1))
		{
			grid_.u()(i - 1, j - 1)   = grid_.u(i - 1, j, EAST) ; 
			grid_.v()(i - 1, j - 1)   = grid_.v(i, j - 1, NORTH) ;
		}

//         if (grid_.isEmpty(i + 1, j + 1))   /////////////////// no values /////////////////////
//             grid_.u()(i, j + 1) = grid_.u()(i, j) + (grid_.dy() / grid_.dx()) * (grid_.v(i + 1, j, WEST) - grid_.v()(i, j)) ;
// 
//         if (grid_.isEmpty(i + 1, j - 1))
//             grid_.u()(i, j - 1) = grid_.u()(i, j) - (grid_.dy() / grid_.dx()) * (grid_.v(i + 1, j - 1, DIAG) - grid_.v(i, j - 1, NORTH)) ;
    }

    // cell(i+1,j) , cell(i,j-1) , cell(i,j+1) are empty ESN
    if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i, j - 1)  && grid_.isEmpty(i, j + 1))
    {
        if (compP)
            grid_.p()(i, j) = 0.0 ;
        grid_.v()(i, j) += gy_ * dt ;
        grid_.v()(i, j - 1) = grid_.v(i, j - 1, NORTH) + gy_ * dt ;
        grid_.u()(i, j) = grid_.u(i - 1, j, EAST) - (grid_.dx() / grid_.dy()) * (grid_.v()(i, j) - grid_.v(i, j - 1, NORTH)) ;

		if (grid_.isEmpty(i + 1, j + 1)) 
		{
			grid_.u()(i, j + 1)     = grid_.u()(i, j) ;
			grid_.v()(i + 1, j)     = grid_.v()(i, j) ;
		}
		if (grid_.isEmpty(i + 1, j - 1)) 
		{
			grid_.u()(i, j - 1)     = grid_.u()(i, j) ; 
			grid_.v()(i + 1, j - 1) = grid_.v(i, j - 1, NORTH) ;
		}

        if (grid_.isEmpty(i - 1, j + 1))
            grid_.u()(i - 1, j + 1) = grid_.u(i - 1, j, EAST) - (grid_.dy() / grid_.dx()) * (grid_.v()(i, j) - grid_.v(i - 1, j, EAST)) ;

        if (grid_.isEmpty(i - 1, j - 1))
            grid_.u()(i - 1, j - 1) = grid_.u(i - 1, j, EAST) + (grid_.dy() / grid_.dx()) * (grid_.v(i, j - 1, NORTH) - grid_.v(i - 1, j - 1, DIAG)) ;
    }
}
//*******************************************************************************************************************

void FluidSimulator::four_empty_neighbour(int i , int j , const real &dt, bool compP)
{
    std::cout << "four_empty_neighbour:" << i << ", " << j << std::endl;
    // According to page 96 of book number 5
    if (grid_.isEmpty(i + 1, j)  && grid_.isEmpty(i - 1, j)  && grid_.isEmpty(i, j - 1)  && grid_.isEmpty(i, j + 1))
    {
        std::cout << "grid_.v()(i, j): " << grid_.v()(i, j) << std::endl;
        std::cout << "grid_.v()(i, j - 1): " << grid_.v()(i, j - 1) << std::endl;
        std::cout << "gy_ * dt: " << gy_ * dt << std::endl;
        if (compP)
            grid_.p()(i, j) = 0.0 ;
        grid_.u()(i, j) += gx_ * dt ;
        grid_.u()(i - 1, j) = grid_.u(i - 1, j, EAST) + gx_ * dt ;
        grid_.v()(i, j) += gy_ * dt ;
        grid_.v()(i, j - 1) = grid_.v(i, j - 1, NORTH) + gy_ * dt ;

        grid_.u()(i, j + 1)     = grid_.u()(i, j) ;
        grid_.u()(i, j - 1)     = grid_.u()(i, j) ;
        grid_.u()(i - 1, j + 1)   = grid_.u(i - 1, j, EAST) ;
        grid_.u()(i - 1, j - 1)   = grid_.u(i - 1, j, EAST) ;

        grid_.v()(i + 1, j)     = grid_.v()(i, j) ;
        grid_.v()(i + 1, j - 1)   = grid_.v(i, j - 1, NORTH) ;
        grid_.v()(i - 1, j)     = grid_.v()(i, j) ;
        grid_.v()(i - 1, j - 1)   = grid_.v(i, j - 1, NORTH) ;
    }
}
