
#include "../src/ParticleTracer.hh"
#include "../src/StaggeredGrid.hh"
#include "../src/Debug.hh"
#include "../src/VTKWriter.hh"

#include <iostream>

int main( )
{
    StaggeredGrid grid(3,3,1,1);
    grid.p().print();
    ParticleTracer tracer(grid);

    tracer.addRectangle(1,1,3,3);

    VTKWriter writer("particletracertest");
    writer.write(grid, &tracer);
   // std::cout << "Copy Test: ";
   // copyTest();
   // std::cout << "OK" << std::endl;

   // std::cout << "Contiguous Memory Test: ";
   // contiguousMemoryTest();
   // std::cout << "OK" << std::endl;

   return 0;
}
