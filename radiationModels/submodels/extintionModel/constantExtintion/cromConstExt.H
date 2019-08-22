/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

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
