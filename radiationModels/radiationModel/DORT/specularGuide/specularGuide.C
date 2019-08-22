/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "specularGuide.H"
#include "fvm.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::specularGuide::specularGuide
(
    const Foam::fvPatch& patch,
    const Foam::List<Foam::vector> mainDirs,
	const List<List< label > > facesNeighbours,
	const List<scalar> faceMagnitudes,
	const List<vector> quadVertex,
	const List<face> quadFaces,
	const bool secondOrder,
	const label nLambda
)
:
    patch_(patch),
    indexes_(patch_.size(),(short)(-1)),
	interpolators_(0),
	storedIntensities_(0)
{
	//stores a specularIndexer for each face direction
	const vectorField n(patch_.nf());
	for(label lambdaI = 0; lambdaI < nLambda; lambdaI++)
	{
		PtrList<scalarField> bandIntensities(mainDirs.size());
		forAll(bandIntensities, dirI)
		{
			bandIntensities.set(dirI, new scalarField(patch_.size(),0.0));
		}
		storedIntensities_.append(bandIntensities);
	}
	forAll(patch_, faceI)
	{
		if(indexes_[faceI]==-1)
		{
			vector currentN = n[faceI];
			indexes_[faceI]=interpolators_.size();
			interpolators_.setSize(interpolators_.size()+1);
			interpolators_.set
			(
				interpolators_.size()-1,
				new specularIndexer
				(
					facesNeighbours,
					mainDirs,
					faceMagnitudes,
					currentN,
					quadVertex,
					quadFaces,
					secondOrder
					
				)
			);
			for(label i=faceI+1;i<patch_.size();i++)
			{
				if(n[i]==currentN) //may be introduced a relaxed criteria
				{
					indexes_[i]=interpolators_.size()-1;
				}
			}
		}
	}
}

// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

void Foam::radiation::specularGuide::updateSpecular
(
	const Foam::PtrList<Foam::scalarField> patchIntensity, const label lambdaI
)
{
	
	forAll(patch_, faceId)
	{
		//generates a scalarField for every face
		scalarField faceIntensities(patchIntensity.size());
		forAll(faceIntensities, dirI)
		{
			faceIntensities[dirI] = patchIntensity[dirI][faceId];
		}
		scalarField reflectedIntensities = 
			interpolators_[indexes_[faceId]].getSpec(faceIntensities);
		forAll(reflectedIntensities, dirI)
		{
			storedIntensities_[lambdaI][dirI][faceId] = reflectedIntensities[dirI];
		}
	}
}

Foam::scalarField Foam::radiation::specularGuide::getSpec
(
	const Foam::label rayId, const Foam::label lambdaId
) const
{
	return storedIntensities_[lambdaId][rayId];
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::specularGuide::~specularGuide()
{}

// ************************************************************************* //
