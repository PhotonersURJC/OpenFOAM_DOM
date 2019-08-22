/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "bandPhaseFunc.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(bandPhaseFunc, 0);

        addToRunTimeSelectionTable
        (
            phaseFuncModel,
            bandPhaseFunc,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::bandPhaseFunc::bandPhaseFunc
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    phaseFuncModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Model")),
    phaseFunctions_(0)
{
	wordList bandNames = coeffsDict_.toc();
	phaseFunctions_.resize(bandNames.size());
	forAll(bandNames, lambdaI)
	{
		dictionary currentBandDict = coeffsDict_.subDict(bandNames[lambdaI]);
		dictionary newBandDict = new dictionary();
		word modelName = currentBandDict.lookupOrDefault<word>("model", "isotropicPhaseFunc");
		newBandDict.add(keyType("phaseFunction"),modelName);
		modelName += "Coeffs";
		newBandDict.add(keyType(modelName),currentBandDict.subDict("modelCoeffs"));
		phaseFunctions_.set
		(
			lambdaI,
			phaseFuncModel::New(newBandDict, mesh_).ptr()
		);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::bandPhaseFunc::~bandPhaseFunc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::bandPhaseFunc::MatrixGen
(const List<vector> directions, const label model)
{
	forAll(phaseFunctions_, lambdaI)
	{
		phaseFunctions_[lambdaI].MatrixGen(directions,model);
	}
}

Foam::scalar Foam::radiation::bandPhaseFunc::sctrFactor
(const label model, const label band, const label directionFrom, const label directionTo)
const
{
	return phaseFunctions_[band].sctrFactor(model, band, directionFrom, directionTo);
}


// ************************************************************************* //
