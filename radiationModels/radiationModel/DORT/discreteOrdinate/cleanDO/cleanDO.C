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

#include "cleanDO.H"
#include "fvm.H"
#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::cleanDO::cleanDO
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
    const List<vector> mainAxis
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
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::cleanDO::~cleanDO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::cleanDO::correct
()
{

    scalar maxResidual = -GREAT;
    forAll(ILambda_, lambdaI)
    {
        const surfaceScalarField Ji(dAve_ & mesh_.Sf());
        fvScalarMatrix IiEq =
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
        );
        IiEq.relax();
        const solverPerformance ILambdaSol = solve
        (
            IiEq,
            "Ii"
        );
        const scalar initialRes =
            ILambdaSol.initialResidual()*omega_/quad_.omegaMax();
		ILambda_[lambdaI].max(0.0);
        maxResidual = max(initialRes, maxResidual);
    }
    return maxResidual;
}

void Foam::radiation::cleanDO::setBandScatter
 (const label lambdaI, const volScalarField value)
{}

void Foam::radiation::cleanDO::changeBand
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
}


// ************************************************************************* //
