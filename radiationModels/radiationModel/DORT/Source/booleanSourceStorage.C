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

#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::booleanSourceStorage::booleanSourceStorage(const label nBands, const label nFaces, const label nDirections)
:
    
    Qlist_(nBands),
    criteria_(nFaces)
{
	forAll(criteria_, critI)
	{
		List<bool> currentCrit(nDirections,false);
		criteria_[critI]=currentCrit;
	}
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::booleanSourceStorage::~booleanSourceStorage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::booleanSourceStorage::setQs(const List<scalar> Qlist)
{
	forAll(Qlist_,lambdaI)
	{
		Qlist_[lambdaI]=Qlist[lambdaI];
	}
}

void Foam::radiation::booleanSourceStorage::setBools
(
	const label face, const label direction, const bool criteria
)
{
	if(criteria)
	{
		criteria_[face][direction]=true;
	}
}

Foam::scalar Foam::radiation::booleanSourceStorage::getValue
(
	const label face, const label direction, const label band
) const
{
	if(criteria_[face][direction])
	{
		return Qlist_[band];
	}
	return 0.0;
}


// ************************************************************************* //
