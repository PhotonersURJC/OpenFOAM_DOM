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

\*---------------------------------------------------------------------------*/

#include "lambertConeRadSource.H"
#include "phaseFuncModel.H"
#include "radExtintionModel.H"
#include "constants.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(lambertConeRadSource, 0);
        addToSourceRunTimeSelectionTables(lambertConeRadSource);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::lambertConeRadSource::initialise()
{
    readIfPresent("Opening", opening_);
	opening_ = (opening_*pi)/180.0;
	scalar directivity = 2.0/(1.0-Foam::cos(opening_*0.5));
	expFactor_ = 0.5*(directivity-2.0);
	correctionFactor_ = (expFactor_+2.0)*0.5/pi;
	opening_=0.0; //avoids quadrature adjust
	sourceValues_.setSize(patchesID_.size());
	//set the direction normal to the first patch in fluxes list
	forAll(mesh_.boundary(), patchI)
	{
		const polyPatch currentPatch(mesh_.boundary()[patchI].patch());
		if(currentPatch.name()==patchesID_[0])
		{
			vector patchDir(0,0,0);
			forAll(currentPatch, faceI)
			{
				patchDir += currentPatch.faceNormals()[faceI];
			}
			direction_=-patchDir/mag(patchDir);
			Info << direction_ << endl;
			break;
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::lambertConeRadSource::lambertConeRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands)
:
    radSource(typeName, mesh, sourceIndex, nBands),
	sourceValues_(0)
{
    initialise();
}


Foam::radiation::lambertConeRadSource::lambertConeRadSource
(
    const dictionary& dict,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
:
    radSource(typeName, dict, mesh, sourceIndex, nBands),
	sourceValues_(0)
{
    initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::lambertConeRadSource::~lambertConeRadSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::lambertConeRadSource::read()
{
    return radSource::read();
}

const Foam::List<Foam::vector> Foam::radiation::lambertConeRadSource::getAxis
(Foam::label nTheta, Foam::label nPhi)
{
    List<vector> mainAxis(3,vector());
    vector Vx(1,0,0); vector Vy(0,1,0); vector Vz(0,0,1);
    if((Vz & direction_) == -1.0)
    { Vx = -Vx; Vy = -Vy; Vz = -Vz; }
    else if ((Vz & direction_) < 1.0)
    {
	vector VCross = Vz ^ direction_;
	VCross = VCross / mag(VCross);
	vector VSum = Vz + direction_;
	VSum = VSum / mag(VSum);
	vector VNormal = VCross ^ VSum;
	Vx = (VCross*VCross.x()) + (VSum*VSum.x()) - (VNormal*VNormal.x());
	Vy = (VCross*VCross.y()) + (VSum*VSum.y()) - (VNormal*VNormal.y());
	Vz = direction_;
    }
	forAll(sourceValues_, patchIndex)
	{
		sourceValues_.set(patchIndex,new thetaSourceStorage(nBands_,nTheta, nPhi));
	}
    mainAxis[0] = Vx; mainAxis[1] = Vy; mainAxis[2] = Vz;
    return mainAxis;
}

Foam::scalar Foam::radiation::lambertConeRadSource::calculate
(
	const word patchID, const label face, const label direction, const label band
) const
{
	if(patchesID_.size()>0)
	{
		forAll(patchesID_, patch)
		{
			if(patchID == patchesID_[patch])
			{
				return sourceValues_[patch].getValue(face, direction, band);
			}
		}
	}
	return 0.0;
}

void Foam::radiation::lambertConeRadSource::generate(const List<vector> directions)
{
	if(patchesID_.size()==0){return;}
    label patchIndex = 0;
    while (patchIndex < patchesID_.size())
    {
        forAll(mesh_.boundary(), patchI)
        {
			const polyPatch currentPatch(mesh_.boundary()[patchI].patch());
            if(currentPatch.name()==patchesID_[patchIndex])
				sourceValues_[patchIndex].setLambertQs(patchesQ_[patchIndex], expFactor_, correctionFactor_);
        }
        patchIndex++;
    }
}




// ************************************************************************* //
