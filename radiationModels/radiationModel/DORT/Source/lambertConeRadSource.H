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
    Foam::radiation::lambertConeRadSource

Description

    lambertian radiation distribution in an hemisphere. Quadrature rotates
	to fit its Z-axis in the source central direction.

SourceFiles
    lambertConeRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceLambertConeRadSource_H
#define radSourceLambertConeRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "thetaSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class lambertConeRadSource Declaration
\*---------------------------------------------------------------------------*/

class lambertConeRadSource
:
    public radSource
{

	PtrList<thetaSourceStorage> sourceValues_;
	scalar expFactor_;
	scalar correctionFactor_;
	

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        lambertConeRadSource(const lambertConeRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const lambertConeRadSource&);


public:

    //- Runtime type information
    TypeName("lambertConeRadSource");


    // Constructors

        //- Construct from components
        lambertConeRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        lambertConeRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~lambertConeRadSource();


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
