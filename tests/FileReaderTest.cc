
#include "../src/FileReader.hh"
#include "../src/Debug.hh"

#include <cmath>


int main( int argc, char** argv )
{

   if ( argc < 2 ) {
       std::cerr << "No config file given" << std::endl;
       return EXIT_FAILURE;
   }

   FileReader myReader;
   
   myReader.registerIntParameter   ("exampleInt" );
   myReader.registerRealParameter  ("exampleReal");
   myReader.registerStringParameter("exampleString" );
   
   CHECK_MSG( myReader.readFile( argv[1] ), "Could not open file " << argv[1] << " which has to be in the current directory." );

   CHECK( myReader.getIntParameter   ("exampleInt")    == 42 );
   CHECK( myReader.getStringParameter("exampleString") == "someStringValue" );
   CHECK( std::abs( myReader.getRealParameter("exampleReal") - 42.4242 ) < 1e-5 );

   myReader.registerIntParameter("aNewInt");
   myReader.setParameter( "aNewInt",    43 ); // add new value ( no registration required )
   myReader.setParameter( "exampleInt", 44 ); // overwrite existing value

   CHECK( myReader.getIntParameter("aNewInt")     == 43 );
   CHECK( myReader.getIntParameter("exampleInt")  == 44 );


   std::cout << "File Reader Test passed successfully" << std::endl;

   return 0;
}
