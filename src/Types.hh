#ifndef TYPES_HH
#define TYPES_HH


// This typedef makes it possible to switch between float and double accuracy
// please do not use "float" or "double" directly in your code, but use real instead
typedef double real;

// Enumeration of boundary conditions
typedef enum { NOSLIP, SLIP, INFLOW, OUTFLOW, PERIODIC } BCTYPE;

// Directions
typedef enum {NORTH, EAST, SOUTH, WEST, DIAG} Direction;

// Flags
typedef unsigned int flag;

const flag OBS = 1 << 0; // OBS = OBSCENTER(old)
const flag FLUID = 1 << 1; // FLUID = FREE(old)

const flag EMPTY = 1 << 2; // EMPTY cell

#endif //TYPES_HH