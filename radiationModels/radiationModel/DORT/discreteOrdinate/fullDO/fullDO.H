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
    Foam::radiation::fullDO

Description
    Discrete Ordinate in an optically heterogenous region
	Full transcription of RTE, extintion coeffcients stored as volFields.

SourceFiles
    fullDO.C

\*---------------------------------------------------------------------------*/

#ifndef fullDO_H
#define fullDO_H

#include "discreteOrdinate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class fullDO Declaration
\*---------------------------------------------------------------------------*/

class fullDO
:
    public discreteOrdinate
{

        //- List of pointers to intensity inScatter fields for given wavelengths
        PtrList<volScalarField> inScatter_;
		
		label actualBand_;

public:

    // Constructors

        //- Construct form components
        fullDO
        (
            const DORT& dort,
	    const quadrature& quad,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label lambda,
            const label rayId,
	    const List<vector> mainAxis
        );


    //- Destructor
    ~fullDO();


    // Member functions

        // Edit

            //- Update radiative intensity on i direction
            scalar correct();

	    void setBandScatter (const label lambdaI, const volScalarField value);
		
		void changeBand (const label lambdaI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
