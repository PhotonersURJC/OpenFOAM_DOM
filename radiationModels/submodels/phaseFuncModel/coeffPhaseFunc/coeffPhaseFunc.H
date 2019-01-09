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
    Foam::radiation::coeffPhaseFunc

Description
    Multi band phase function, same model in all bands, different coeffs.

SourceFiles
    coeffPhaseFunc.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCoeffPhaseFunc_H
#define radiationCoeffPhaseFunc_H

#include "phaseFuncModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                        Class coeffPhaseFunc Declaration
\*---------------------------------------------------------------------------*/

class coeffPhaseFunc
:
    public phaseFuncModel
{
    // Private data

        //- PhaseFunc model dictionary
        dictionary coeffsDict_;
		
		//- List of band phase functions
		PtrList<phaseFuncModel> phaseFunctions_;

public:

    //- Runtime type information
    TypeName("coeffPhaseFunction");


    // Constructors

        //- Construct from components
        coeffPhaseFunc(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~coeffPhaseFunc();

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
