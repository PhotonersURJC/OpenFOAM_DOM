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

#include "thetaSourceStorage.H"
#include "constants.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::thetaSourceStorage::thetaSourceStorage(const label nBands, const label nTheta, const label nPhi)
:
    
    Qlist_(nBands),
    nTheta_(nTheta/2),
	nPhi_(nPhi)
{
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::thetaSourceStorage::~thetaSourceStorage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::thetaSourceStorage::setLambertQs(const List<scalar> Qlist, const scalar expFactor, const scalar correctionFactor)
{
	List<scalar> thetaFactors(nTheta_);
	scalar deltaTheta = (Foam::constant::mathematical::pi*0.5)/nTheta_;
	Info << "nTheta: " << nTheta_ << endl;
	Info << "deltaTheta: " << deltaTheta << endl;
	forAll(thetaFactors, thetaI)
	{
		thetaFactors[thetaI] = correctionFactor*Foam::pow(Foam::cos((thetaI+0.5)*deltaTheta),expFactor);
		Info << "thetaFactors[" << thetaI << "]: " << thetaFactors[thetaI] << endl;
	}
	forAll(Qlist_,lambdaI)
	{
		List<scalar> bandQs(nTheta_);
		forAll(thetaFactors, thetaI)
		{
			bandQs[thetaI]=Qlist[lambdaI]*thetaFactors[thetaI];
		}
		Qlist_[lambdaI]=bandQs;
	}
}

Foam::scalar Foam::radiation::thetaSourceStorage::getValue
(
	const label face, const label direction, const label band
) const
{
	label thetaI = direction/(4*nPhi_);
	if(thetaI<nTheta_)
		return Qlist_[band][thetaI];
	return 0.0;
}


// ************************************************************************* //
