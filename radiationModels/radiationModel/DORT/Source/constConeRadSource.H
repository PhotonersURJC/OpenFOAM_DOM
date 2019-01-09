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
    Foam::radiation::constConeRadSource

Description

    Cone-shaped emmision source. Quadrature rotates to fit its Z-axis in the
	cone central direction, and sets the nearest cut between ordinates in the
	cone opening perimeter. Theta increments in each side of the cut are
	corrected to be equally spaced in their respective regions.

SourceFiles
    constConeRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceconstConeRadSource_H
#define radSourceconstConeRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class constConeRadSource Declaration
\*---------------------------------------------------------------------------*/

class constConeRadSource
:
    public radSource
{

	PtrList<booleanSourceStorage> sourceValues_;
    List<vector> patchNormals_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        constConeRadSource(const constConeRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const constConeRadSource&);


public:

    //- Runtime type information
    TypeName("constConeRadSource");


    // Constructors

        //- Construct from components
        constConeRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        constConeRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~constConeRadSource();


    // Member functions

        // Edit

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
