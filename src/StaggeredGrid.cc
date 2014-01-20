
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

    Array pp( (int) (xxSize / ddx) + 2, (int) (yySize / ddy) + 2 );
    Array rhss( (int) (xxSize / ddx), (int) (yySize / ddy) );
    Array uu( (int) (xxSize / ddx) + 1, (int) (yySize / ddy) + 2 );
    Array vv( (int) (xxSize / ddx) + 2, (int) (yySize / ddy) + 1 );
    Array ff( (int) (xxSize / ddx) + 1, (int) (yySize / ddy) );
    Array gg( (int) (xxSize / ddx), (int) (yySize / ddy) + 1 );
    FlagArray ob( (int) (xxSize / ddx) + 2, (int) (yySize / ddy) + 2 );

    p_ = pp;
    rhs_ = rhss;
    u_ = uu;
    v_ = vv;
    f_ = ff;
    g_ = gg;

    //    Flag array:
    //    obstacle cell (center): ob(i,j) =  1
    //    obstacle cell (west):   ob(i,j) =  2
    //    obstacle cell (east):   ob(i,j) =  4
    //    obstacle cell (north):  ob(i,j) =  8
    //    obstacle cell (south):  ob(i,j) = 16
    ob.fill(FREE);
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
    int imax = configuration.getIntParameter("imax");
    CHECK_MSG( (imax >= 0), "wrong input for imax: " << imax);
    int jmax = configuration.getIntParameter("jmax");
    CHECK_MSG( (jmax >= 0), "wrong input for jmax: " << jmax);

    dx_ = xl / imax;
    dy_ = yl / jmax;
    xSize_ = xl;
    ySize_ = yl;
    Array pp( imax + 2, jmax + 2 );
    Array rhss( imax, jmax );
    Array uu( imax + 1, jmax + 2 );
    Array vv( imax + 2, jmax + 1 );
    Array ff( imax + 1, jmax  );
    Array gg( imax, jmax + 1 );
    FlagArray ob( imax + 2, jmax + 2 );

    p_ = pp;
    rhs_ = rhss;
    u_ = uu;
    v_ = vv;
    f_ = ff;
    g_ = gg;

    //    Flag array:
    //    obstacle cell (center): ob(i,j) =  1
    //    obstacle cell (west):   ob(i,j) =  2
    //    obstacle cell (east):   ob(i,j) =  4
    //    obstacle cell (north):  ob(i,j) =  8
    //    obstacle cell (south):  ob(i,j) = 16
    ob.fill(FREE);
    obs_ = ob;

}

void StaggeredGrid::setCellToObstacle(int x, int y)
{
    obs_(x, y) = obs_(x, y) | OBSCENTER; // set cell to obstacle
    if ( x < u_.getSize(0) && y < u_.getSize(1) )
        u_(x, y) = 0;
    if ( x < v_.getSize(0) && y < v_.getSize(1) )
        v_(x, y) = 0;

    // configurate neighbors
    if ( (x - 1) >= 0 ) // east
        obs_(x - 1, y) = obs_(x - 1, y) | OBSWEST;
    if ( (x + 1) < obs_.getSize(0) ) // west
        obs_(x + 1, y) = obs_(x + 1, y) | OBSEAST;
    if ( (y - 1) >= 0 ) // south
        obs_(x, y - 1) = obs_(x, y - 1) | OBSNORTH;
    if ( (y + 1) < obs_.getSize(1) ) // north
        obs_(x, y + 1) = obs_(x, y + 1) | OBSSOUTH;
}

void StaggeredGrid::createRectangle(real x1, real y1, real x2, real y2)
{
    // define corners in grid
    int startX = (int) (x1 / dx_);
    int startY = (int) (y1 / dy_);
    int endX = (int) (x2 / dx_);
    int endY = (int) (y2 / dy_);

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

            if ( point <= r ) // check if point is in the circle
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
    ob.fill(FREE);

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
