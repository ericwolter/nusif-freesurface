
#include "StaggeredGrid.hh"

StaggeredGrid::StaggeredGrid( int xxSize, int yySize, real ddx, real ddy )
{
    CHECK_MSG( (xxSize >= 0), "wrong input for xSize: " << xxSize);
    CHECK_MSG( (yySize >= 0), "wrong input for ySize: " << yySize);
    CHECK_MSG( (ddx >= 0), "wrong input for dx: " << ddx);
    CHECK_MSG( (ddy >= 0), "wrong input for dy: " << ddy);
    PROG("construct grid manually");
    // initilize grid values
    dx_ = ddx;
    dy_ = ddy;
    xSize_ = xxSize;
    ySize_ = yySize;
    imax_ = (int) (xSize_ / dx_);
    jmax_ = (int) (ySize_ / dy_);

    Array pp( imax_ + 2, jmax_ + 2 );
    Array rhss( imax_, jmax_ );
    Array uu( imax_ + 1, jmax_ + 2 );
    Array vv( imax_ + 2, jmax_ + 1 );
    Array ff( imax_ + 1, jmax_ );
    Array gg( imax_, jmax_ + 1 );
    FlagArray ob( imax_ + 2, jmax_ + 2 );

    p_ = pp;
    rhs_ = rhss;
    u_ = uu;
    v_ = vv;
    f_ = ff;
    g_ = gg;

    //    Flag array:
    //    obstacle cell (center): ob(i,j) =  1
    //    fluid cell (center):    ob(i,j) =  2

    //    empty cell (center):    ob(i,j) =  4
    //    empty cell (north):     ob(i,j) =  8
    //    empty cell (south):     ob(i,j) = 16
    //    empty cell (west):      ob(i,j) = 32
    //    empty cell (east):      ob(i,j) = 64
    ob.fill(FLUID);
    obs_ = ob;

}

StaggeredGrid::StaggeredGrid( const FileReader &configuration )
{
    PROG("construct grid with file");
    // initilize grid values

    real xl = configuration.getRealParameter("xlength");
    CHECK_MSG( (xl >= 0), "wrong input for xlength: " << xl);
    real yl = configuration.getRealParameter("ylength");
    CHECK_MSG( (yl >= 0), "wrong input for ylength: " << yl);
    int iimax = configuration.getIntParameter("imax");
    CHECK_MSG( (iimax >= 0), "wrong input for imax: " << iimax);
    int jjmax = configuration.getIntParameter("jmax");
    CHECK_MSG( (jjmax >= 0), "wrong input for jmax: " << jjmax);

    dx_ = xl / iimax;
    dy_ = yl / jjmax;
    xSize_ = xl;
    ySize_ = yl;
    imax_ = iimax;
    jmax_ = jjmax;
    Array pp( imax_ + 2, jmax_ + 2 );
    Array rhss( imax_, jmax_ );
    Array uu( imax_ + 1, jmax_ + 2 );
    Array vv( imax_ + 2, jmax_ + 1 );
    Array ff( imax_ + 1, jmax_  );
    Array gg( imax_, jmax_ + 1 );
    FlagArray ob( imax_ + 2, jmax_ + 2 );

    p_ = pp;
    rhs_ = rhss;
    u_ = uu;
    v_ = vv;
    f_ = ff;
    g_ = gg;

    //    Flag array:
    //    obstacle cell (center): ob(i,j) =  1
    //    fluid cell (center):    ob(i,j) =  2

    //    empty cell (center):    ob(i,j) =  4
    //    empty cell (north):     ob(i,j) =  8
    //    empty cell (south):     ob(i,j) = 16
    //    empty cell (west):      ob(i,j) = 32
    //    empty cell (east):      ob(i,j) = 64
    ob.fill(FLUID);
    obs_ = ob;

}

void StaggeredGrid::createRectangle(real x1, real y1, real x2, real y2)
{
    // define corners in grid
    int startX = std::min( (int) (x1/dx_), (int) (x2/dx_) );
    int endX = std::max( (int) (x1/dx_), (int) (x2/dx_) );
    int startY = std::min( (int) (y1/dy_), (int) (y2/dy_) );
    int endY = std::max( (int) (y1/dy_), (int) (y2/dy_) );

    for ( int i = startX; i <= endX; ++i )
    {
        for ( int j = startY; j <= endY; ++j )
        {
            //PROG("in for");
            setCellToObstacle(i, j); // set cell to obstacle
        }
    }

}

void StaggeredGrid::createCircle(real x, real y, real r)
{
    real point = 0.0;
    for ( int i = 0; i < obs_.getSize(0); ++i )
    {
        for ( int j = 0; j < obs_.getSize(1); ++j )
        {

            // compute vector to the actual point
            point = (x - i * dx_) * (x - i * dx_) + (y - j * dy_) * (y - j * dy_);

            if ( point <= r*r ) // check if point is in the circle
                setCellToObstacle(i, j); // set cell to obstacle
        }
    }

}

void StaggeredGrid::createPng( const std::string &pngFilename )
{
    // create png file
    int n = obs_.getSize(0);
    int m = obs_.getSize(1);
    GrayScaleImage image(pngFilename, n, m);

    // set size
    image = image.getResizedImage( n, m );
    for ( int i = 0; i < n; ++i )
    {
        for ( int j = 0; j < m; ++j )
        {
            if ( isFluid(i, j) )  // set cell to fluid
            {
                image(i, j) = std::numeric_limits<unsigned char>::max();
            }
            else     // set cell to obstacle
            {
                image(i, j) = 0;
            }
        }
    }

    // save image
    image.save( pngFilename );

}

void StaggeredGrid::readPng( const std::string &pngFilename )
{
    // load png file
    GrayScaleImage image( pngFilename );

    // get size
    int n = image.width();
    int m = image.height();

    // create flag array
    FlagArray ob( n, m );
    ob.fill(FLUID);

    // assign to obstacle array
    obs_ = ob;
    for ( int i = 0; i < n; ++i )
    {
        for ( int j = 0; j < m; ++j )
        {
            if ( image.getElement(i, j) == 0 )  // set cell to obstacle
            {
                setCellToObstacle(i, j);
            }
        }
    }

}


real StaggeredGrid::u_inter ( real x , real y )
{
    int i = (int) (x / dx_) + 1                 ;
    int j = (int) (( y + 0.5 * dy_ ) / dy_) + 1   ;

    real x1 , x2 , y1 , y2 ;

    x1 = (real) (i - 1) * dx_       ;
    x2 = (real) i * dx_           ;
    y1 = (real) ((j - 1) - 0.5) * dy_ ;
    y2 = (real) (j - 1) * dy_       ;

    real u1 = u_( i - 1 , j - 1 ) ;
    real u2 = u_( i   , j - 1 ) ;
    real u3 = u_( i - 1 , j   ) ;
    real u4 = u_( i   , j   ) ;

    real u_interpolate ;

    real coeeficient = (real) ( 1.0 / (dx_ * dy_) ) ;

    u_interpolate = coeeficient * (   (real) ( x2 - x  ) * ( y2 - y  ) * u1
                                      + (real) ( x  - x1 ) * ( y2 - y  ) * u2
                                      + (real) ( x2 - x  ) * ( y  - y1 ) * u3
                                      + (real) ( x  - x1 ) * ( y  - y1 ) * u4
                                  ) ;

    return u_interpolate ;

}

real StaggeredGrid::v_inter ( real x , real y )
{
    int i = (int) (( x + 0.5 * dx_ ) / dx_) + 1 ;
    int j = (int) ( y  / dy_) + 1               ;

    real x1 , x2 , y1 , y2 ;

    x1 = (real) ((i - 1) - 0.5) * dx_ ;
    x2 = (real) (i - 0.5)    * dx_ ;
    y1 = (real) (j - 1) * dy_       ;
    y2 = (real)     j * dy_       ;

    real v1 = v_( i - 1 , j - 1 ) ;
    real v2 = v_( i   , j - 1 ) ;
    real v3 = v_( i - 1 , j   ) ;
    real v4 = v_( i   , j   ) ;

    real v_interpolate ;

    real coeeficient = (real) ( 1.0 / (dx_ * dy_) ) ;

    v_interpolate = coeeficient * (   (real) ( x2 - x  ) * ( y2 - y  ) * v1
                                      + (real) ( x  - x1 ) * ( y2 - y  ) * v2
                                      + (real) ( x2 - x  ) * ( y  - y1 ) * v3
                                      + (real) ( x  - x1 ) * ( y  - y1 ) * v4
                                  ) ;

    return v_interpolate;
}

void StaggeredGrid::setCellToFluid(int x, int y)
{
    obs_(x,y) = FLUID; // set cell to fluid
}

void StaggeredGrid::setCellToEmpty(int x, int y)
{
    obs_(x, y) = EMPTY; // set cell to empty
    if ( x < u_.getSize(0) && y < u_.getSize(1) )
        u_(x, y) = 0;
    if ( x < v_.getSize(0) && y < v_.getSize(1) )
        v_(x, y) = 0;
}

void StaggeredGrid::setCellToObstacle(int x, int y)
{
    obs_(x, y) = OBS; // set cell to obstacle
    if ( x < u_.getSize(0) && y < u_.getSize(1) )
        u_(x, y) = 0;
    if ( x < v_.getSize(0) && y < v_.getSize(1) )
        v_(x, y) = 0;

//     // configurate neighbors
//     if ( (x - 1) >= 0 ) // east
//         obs_(x - 1, y) = obs_(x - 1, y) | OBSWEST;
//     if ( (x + 1) < obs_.getSize(0) ) // west
//         obs_(x + 1, y) = obs_(x + 1, y) | OBSEAST;
//     if ( (y - 1) >= 0 ) // south
//         obs_(x, y - 1) = obs_(x, y - 1) | OBSNORTH;
//     if ( (y + 1) < obs_.getSize(1) ) // north
//         obs_(x, y + 1) = obs_(x, y + 1) | OBSSOUTH;
}

void StaggeredGrid::markCells()
{
    for (int i = 1; i <= imax_; ++i)
    {
        for (int j = 1; j <= jmax_; ++j)
        {
            if (!isFluid(i, j)) continue;

            setCellToEmpty(i, j);
        }
    }

    for (std::vector<Particle>::iterator it = particles_.begin() ; it != particles_.end(); ++it)
    {
        Particle prtcl = *it;

        int i = prtcl.getCellX(dx_);
        int j = prtcl.getCellY(dy_);

        setCellToFluid(i, j);
    }
}

