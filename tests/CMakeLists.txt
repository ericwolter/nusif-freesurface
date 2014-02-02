
include_directories( . )

add_executable( tracer_test ../src/FileReader.cc ../src/Array.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/ParticleTracer.cc ../src/SORSolver.cc ../src/VTKWriter.cc ../src/FluidSimulator.cc  ../src/Particle.cc ../src/GrayScaleImage.cc ../src/lodepng.cc  ParticleTracerTests.cc )

add_executable( surface_test ../src/FileReader.cc ../src/Array.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/ParticleTracer.cc ../src/SORSolver.cc ../src/VTKWriter.cc ../src/FluidSimulator.cc  ../src/Particle.cc ../src/GrayScaleImage.cc ../src/lodepng.cc  FreeSurfaceTests.cc )

