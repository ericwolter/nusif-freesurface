
include_directories( . )

add_executable( sor_test ../src/FileReader.cc ../src/Array.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/VTKWriter.cc ../src/SORSolver.cc ../src/GrayScaleImage.cc ../src/lodepng.cc SORTest.cc )

add_executable( array_test      ../src/Array.cc      ../src/Debug.cc ArrayTest.cc      )
add_executable( filereader_test ../src/FileReader.cc ../src/Debug.cc FileReaderTest.cc )

add_executable( fg_test ../src/FileReader.cc ../src/Array.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc ../src/VTKWriter.cc ../src/FluidSimulator.cc ../src/GrayScaleImage.cc ../src/lodepng.cc FGTest.cc )

add_executable( dcavity_test ../src/FileReader.cc ../src/Array.cc ../src/Debug.cc ../src/StaggeredGrid.cc ../src/SORSolver.cc ../src/VTKWriter.cc ../src/FluidSimulator.cc ../src/GrayScaleImage.cc ../src/lodepng.cc  DcavityTest.cc )

