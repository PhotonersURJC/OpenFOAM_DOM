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
    Foam::radiation::isotropicRadSource

Description

    Isotropic emmision source. Light is equally distributed in all discrete
	ordinates with their central direction leaving the face. Quadrature keeps
	the mesh coordinates system.

SourceFiles
    isotropicRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceisotropicRadSource_H
#define radSourceisotropicRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class isotropicRadSource Declaration
\*---------------------------------------------------------------------------*/

class isotropicRadSource
:
    public radSource
{
	// Private data
		PtrList<booleanSourceStorage> sourceValues_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        isotropicRadSource(const isotropicRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const isotropicRadSource&);


public:

    //- Runtime type information
    TypeName("isotropicRadSource");


    // Constructors

        //- Construct from components
        isotropicRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        isotropicRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~isotropicRadSource();


    // Member functions

	    //- Return main axis to draw the new quadrature
	    const List<vector> getAxis(label nTheta, label nPhi);
		
            //- Store values of emission
            virtual void generate(const List<vector> directions);

            //- Read radiation properties dictionary
            bool read();

        // Access
		
			//- Acces to specific stored value
            scalar calculate(const word patchID, const label face, const label direction, const label band) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
