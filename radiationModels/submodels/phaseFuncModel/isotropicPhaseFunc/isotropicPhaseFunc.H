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
    Foam::radiation::isotropicPhaseFunc

Description
    Isotropic phase function

SourceFiles
    isotropicPhaseFunc.C

\*---------------------------------------------------------------------------*/

#ifndef radiationIsotropicPhaseFunc_H
#define radiationIsotropicPhaseFunc_H

#include "phaseFuncModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class isotropicPhaseFunc Declaration
\*---------------------------------------------------------------------------*/

class isotropicPhaseFunc
:
    public phaseFuncModel
{
    // Private data

        //- Scattering isotropic factor
		scalar factor_;

public:

    //- Runtime type information
    TypeName("isotropicPhaseFunction");


    // Constructors

        //- Construct from components
        isotropicPhaseFunc(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~isotropicPhaseFunc();

	    //Scattering Matrix Generation
	    void MatrixGen(const List<vector> directions, const label model);
		
		scalar sctrFactor(const label model, const label band, const label directionFrom, const label directionTo) const;
				
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
