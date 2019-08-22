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

#include "thermoAbsorptionDO.H"
#include "fvm.H"
#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::thermoAbsorptionDO::thermoAbsorptionDO
(
    const DORT& dort,
    const quadrature& quad,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const label rayId,
    const List<vector> mainAxis,
    const scalar k
)
:
    discreteOrdinate
    (
        dort,
        quad,
        mesh,
        phi,
        theta,
        deltaPhi,
        deltaTheta,
        nLambda,
        rayId,
        mainAxis
    ),
    k_("k", dimless/dimLength, k)
{
    k_ = k_*omega_; //To avoid repeat it in every iteration
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::thermoAbsorptionDO::~thermoAbsorptionDO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::thermoAbsorptionDO::correct
()
{
    scalar maxResidual = -GREAT;
    forAll(ILambda_, lambdaI)
    {
        const surfaceScalarField Ji(dAve_ & mesh_.Sf());
		const volScalarField emission = dort_.Elambda(lambdaI);
		if(nLambda_>1)
		{
			k_.value()=max(dort_.kappaLambda(lambdaI)).value()*omega_;
		}
        fvScalarMatrix IiEq =
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
          + fvm::Sp(k_, ILambda_[lambdaI]) == emission
		);
        IiEq.relax();
        const solverPerformance ILambdaSol = solve
        (
            IiEq,
            mesh_.solver("Ii")
        );
        const scalar initialRes =
            ILambdaSol.initialResidual()*omega_/quad_.omegaMax();
		ILambda_[lambdaI].max(0.0);
        maxResidual = max(initialRes, maxResidual);
    }

    return maxResidual;
}

void Foam::radiation::thermoAbsorptionDO::setBandScatter
 (const label lambdaI, const volScalarField value)
 {}

void Foam::radiation::thermoAbsorptionDO::changeBand
 (const label lambdaI)
{
	IOobject IHeader
    (
        "IBand_" +  name(lambdaI),
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
	if (IHeader.good())
	{
		ILambda_.set
		(
			0,
			new volScalarField(IHeader, mesh_)
		);
	}
	k_=max(dort_.kappaLambda(lambdaI)).value()*omega_;
}


// ************************************************************************* //
