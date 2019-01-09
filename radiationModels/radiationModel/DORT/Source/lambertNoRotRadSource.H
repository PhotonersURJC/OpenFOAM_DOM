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
    Foam::radiation::lambertNoRotRadSource

Description

    lambertian radiation distribution in an hemisphere.
	No adjust in quadrature

SourceFiles
    lambertNoRotRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceLambertNoRotRadSource_H
#define radSourceLambertNoRotRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "thetaSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class lambertNoRotRadSource Declaration
\*---------------------------------------------------------------------------*/

class lambertNoRotRadSource
:
    public radSource
{
	scalar expFactor_;
	scalar correctionFactor_;
	List<List<List<scalar> > > sourceValues_;
	

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        lambertNoRotRadSource(const lambertNoRotRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const lambertNoRotRadSource&);


public:

    //- Runtime type information
    TypeName("lambertNoRotRadSource");


    // Constructors

        //- Construct from components
        lambertNoRotRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        lambertNoRotRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~lambertNoRotRadSource();


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
