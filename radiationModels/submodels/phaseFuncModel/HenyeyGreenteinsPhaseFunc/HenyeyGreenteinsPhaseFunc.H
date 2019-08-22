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
