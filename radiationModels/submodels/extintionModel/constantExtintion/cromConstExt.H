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
    Foam::radiation::cromConstExtintion

Description
    Homogeneous radiation extintion coefficient for multi-band
	
	\verbatim
		extintionModel cromConstExt;
		cromConstExtCoeffs
		{
		   absorptivity
		   {
			   band1	absorptivity    [ m^-1 ]         0.2;
			   band2	absorptivity    [ m^-1 ]         0.5;
		   }
		   scatter
		   {
			   band1	scatter    [ m^-1 ]         0.3;
			   band2	scatter    [ m^-1 ]         0.4;
		   }
		   bandLimits
		   {
			   band1 (300 400);
			   band2 (400 500);
		   }
		}
	\endverbatim

SourceFiles
    cromConstExt.C

\*---------------------------------------------------------------------------*/

#ifndef radiationCromConstExt_H
#define radiationCromConstExt_H

#include "radExtintionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                  Class cromConstExt Declaration
\*---------------------------------------------------------------------------*/

class cromConstExt
:
    public radExtintionModel
{
    // Private data

        //- Absorption model dictionary
        dictionary coeffsDict_;

        //- Band absorption coefficient / [1/m]
        PtrList<dimensionedScalar> kappaLambda_;
		
		//- Band scattering coefficient / [1/m]
        PtrList<dimensionedScalar> sigmaLambda_;
		
		//- Number of wavelength bands
		label nBands_;
		
		//- Wavelength limits of bands
		List< Vector2D < scalar > > bandLimits_;

public:

    //- Runtime type information
    TypeName("cromConstExt");


    // Constructors

        //- Construct from components
        cromConstExt(const dictionary& dict, const fvMesh& mesh);


    //- Destructor
    virtual ~cromConstExt();


    // Member Functions

        // Access

            // Coefficients

                //- Absorption coefficient for continuous phase
                tmp<volScalarField> kappaCont(const label bandI) const;

                //- Scattering coefficient for continuous phase
                tmp<volScalarField> sigmaCont(const label bandI) const;
				
				virtual List< Vector2D < scalar > > getBandsLimits() const;
		

    // Member Functions

	        //- Correct absorption coefficients
            virtual void kappaCorrect
            (
                volScalarField& kappa,
                PtrList<volScalarField>& kappaj
            ) const;

			//- Correct scattering coefficients
            virtual void sigmaCorrect
            (
                volScalarField& sigma,
                PtrList<volScalarField>& sigmaj
            ) const;
	
        inline bool isGrey() const
        {
            return false;
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
