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
    Foam::radiation::constantExtintion

Description
    Homogeneous radiation extintion coefficient for grey simulations
	
	\verbatim
		extintionModel greyConstExt;
		greyConstExtCoeffs
		{
		   absorptivity    absorptivity    [ m^-1 ]         1.0;
		   scatter         scatter         [ m^-1 ]         0.2;
		}
	\endverbatim

SourceFiles
    constantExtintion.C

\*---------------------------------------------------------------------------*/

#ifndef radiationGreyConstExt_H
#define radiationGreyConstExt_H

#include "radExtintionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class greyConstExt Declaration
\*---------------------------------------------------------------------------*/

class greyConstExt
:
    public radExtintionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Absorption coefficient / [1/m]
        dimensionedScalar kappa_;
		
		//- Absorption coefficient / [1/m]
        dimensionedScalar sigma_;

public:

    //- Runtime type information
    TypeName("greyConstExt");


    // Constructors

        //- Construct from components
        greyConstExt(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~greyConstExt();


    // Member Functions

        // Access

            // Absorption coefficient

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> kappaCont(const label bandI = 0) const;

                //- Scattering coefficient for continuous phase
                tmp<volScalarField> sigmaCont(const label bandI = 0) const;
		

    // Member Functions

        inline bool isGrey() const
        {
            return true;
        }

    inline bool isHomogeneous() const
    { return true;}
    inline bool isConstant() const
    { return true;}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
