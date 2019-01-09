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
    Foam::radiation::HenyeyGreenteinsPhaseFunc

Description
    Henyey-Greenteins phase function

SourceFiles
    HenyeyGreenteinsPhaseFunc.C

\*---------------------------------------------------------------------------*/

#ifndef radiationHenyeyGreenteinsPhaseFunc_H
#define radiationHenyeyGreenteinsPhaseFunc_H

#include "phaseFuncModel.H"
#include "sctrMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
	
	class sctrMatrix;

/*---------------------------------------------------------------------------*\
                  Class HenyeyGreenteinsPhaseFunc Declaration
\*---------------------------------------------------------------------------*/

class HenyeyGreenteinsPhaseFunc
:
    public phaseFuncModel
{
    // Private data


        //- PhaseFunc model dictionary
        dictionary coeffsDict_;

		dimensionedScalar g_;
		
		//3D matrix to store multiband scattering factors;
		List< List < List < scalar > > > sctrMatrix_;
		
		label matrixSize_;

public:

    //- Runtime type information
    TypeName("HenyeyGreenteinsPhaseFunction");


    // Constructors

        //- Construct from components
        HenyeyGreenteinsPhaseFunc(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~HenyeyGreenteinsPhaseFunc();

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
