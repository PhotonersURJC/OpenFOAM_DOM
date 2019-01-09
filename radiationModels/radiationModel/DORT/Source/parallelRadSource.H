/*---------------------------------------------------------------------------*\
             Discrete Ordinate Method Radiation Model for OpenFOAM.


Code corresponding to the article entitled:


"Improved Discrete Ordinate Method for accurate simulation radiation transport
                  using solar and LED light sources"


by:


José Moreno, Cintia Casado, Javier Marugán

Department of Chemical and Environmental Technology,

ESCET, Universidad Rey Juan Carlos,

C/Tulipán s/n, 28933 Móstoles (Madrid), Spain

Tel. +34 91 664 7466; E-mail: javier.marugan@urjc.es

Class
    Foam::radiation::parallelRadSource

Description

    Parallel emmision source. Light is emmited in a single discrete ordinate.
	Quadrature is rotated to place a discrete ordinate near the quadrature
	equator in the light direction.
	
	Discrete ordinates near the equator are the biggest, and have a squared
	shape, been considered as "better" ordinates.

SourceFiles
    parallelRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceparallelRadSource_H
#define radSourceparallelRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class parallelRadSource Declaration
\*---------------------------------------------------------------------------*/

class parallelRadSource
:
    public radSource
{
	// Private data
	PtrList<booleanSourceStorage> sourceValues_;
    List<vector> patchNormals_;
    label quadIndex_;
    scalar omega_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        parallelRadSource(const parallelRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const parallelRadSource&);


public:

    //- Runtime type information
    TypeName("parallelRadSource");


    // Constructors

        //- Construct from components
        parallelRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        parallelRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~parallelRadSource();


    // Member functions
	
	    //- Return main axis to draw the new quadrature
	    const List<vector> getAxis(label nTheta, label nPhi);

		    //- Store values of emission
            virtual void generate(const List<vector> directions);

            //- Read radiation properties dictionary
            bool read();

        // Access

			//- Acces to specific stored value
            virtual scalar calculate(const word patchID, const label face, const label direction, const label band) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
