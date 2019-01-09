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
    Foam::radiation::sctrMatrix

Description
    Struct to store scattering factors between a pair of directions.

SourceFiles
    sctrMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef sctrMatrix_H
#define sctrMatrix_H

#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class sctrMatrix Declaration
\*---------------------------------------------------------------------------*/

class sctrMatrix
{

        //- List of specular index in patch faces
        List<scalar> sctrFactors_;
		
		label directions_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const sctrMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const sctrMatrix&);


public:

    // Constructors

        //- Construct from components
        sctrMatrix
        (
            const label directions
        );
		


    //- Destructor
    ~sctrMatrix();


    // Member functions

		//- Return non-const access to sctrFactor from dir i to dir j
			scalar& Item(const label i, const label j);
			
		//- Return const access to sctrFactor from dir i to dir j
		scalar constItem(const label i, const label j) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
