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

#include "sctrMatrix.H"
#include "fvm.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::sctrMatrix::sctrMatrix
(
    const Foam::label directions
)
:
    sctrFactors_(directions*directions),
	directions_(directions)
{

}

Foam::scalar& Foam::radiation::sctrMatrix::Item
(const label i, const label j)
{
	return sctrFactors_[i*directions_+j];
}

Foam::scalar Foam::radiation::sctrMatrix::constItem
(const label i, const label j) const
{
	return sctrFactors_[i*directions_+j];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::sctrMatrix::~sctrMatrix()
{}

// ************************************************************************* //
