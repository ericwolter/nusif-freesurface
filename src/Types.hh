#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;

// Enumeration of boundary conditions
typedef enum { NOSLIP, SLIP, INFLOW, OUTFLOW, PERIODIC } BCTYPE;

// Directions
typedef enum {NORTH, EAST, SOUTH, WEST} Direction;

// Flags
typedef unsigned char flag;

const flag OBSCENTER = 1 << 0;
const flag OBSWEST   = 1 << 1;
const flag OBSEAST   = 1 << 2;
const flag OBSNORTH  = 1 << 3;
const flag OBSSOUTH  = 1 << 4; 
const flag FREE      = 1 << 5;

#endif //TYPES_HH
