/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::radiation::thermoAbsorptionDO

Description
    Discrete Ordinate in a region with homogeneous optical properties
	Scattering matrix includes scattering coefficent, to optimize iterations

SourceFiles
    thermoAbsorptionDO.C

\*---------------------------------------------------------------------------*/

#ifndef thermoAbsorptionDO_H
#define thermoAbsorptionDO_H

#include "discreteOrdinate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class absorptionDO Declaration
\*---------------------------------------------------------------------------*/

class thermoAbsorptionDO
:
    public discreteOrdinate
{

	//- absorption coefficient
    dimensionedScalar k_;
	
	//- Black body reference
    const blackBodyEmissionRev& blackBody_;

public:

    // Constructors

        //- Construct form components
        thermoAbsorptionDO
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
	    const List<vector> mainAxis,
            const scalar k,
			const blackBodyEmissionRev& blackBody
        );


    //- Destructor
    ~thermoAbsorptionDO();


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
