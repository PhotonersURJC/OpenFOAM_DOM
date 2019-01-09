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
    Foam::radiation::cleanDO

Description
    Discrete ordinate in transparent region.
	RTE is reduced to Divergence(Radiative Intensity) = 0

SourceFiles
    cleanDO.C

\*---------------------------------------------------------------------------*/

#ifndef cleanDO_H
#define cleanDO_H

#include "discreteOrdinate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class cleanDO Declaration
\*---------------------------------------------------------------------------*/

class cleanDO
:
    public discreteOrdinate
{

public:

    // Constructors

        //- Construct form components
        cleanDO
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
    ~cleanDO();


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
