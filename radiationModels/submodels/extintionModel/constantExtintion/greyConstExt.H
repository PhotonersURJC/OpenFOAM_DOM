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
