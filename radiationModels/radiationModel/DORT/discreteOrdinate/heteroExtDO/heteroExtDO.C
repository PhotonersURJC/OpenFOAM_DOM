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

#include "heteroExtDO.H"
#include "fvm.H"
#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::heteroExtDO::heteroExtDO
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
    ),
    inScatter_(nLambda),
	actualBand_(0)
{
    forAll(ILambda_, lambdaI)
    {
	inScatter_.set
	(
	    lambdaI,
	    new volScalarField
	    (
		IOobject
		(
		    "inScatter",
		    mesh_.time().timeName(),
		    mesh_,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("inScatter", dimMass/pow3(dimTime), 0.0)
	    )
	);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::heteroExtDO::~heteroExtDO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::heteroExtDO::correct
()
{
    scalar maxResidual = -GREAT;

    forAll(ILambda_, lambdaI)
    {
		label actualLambda=lambdaI;
		if(actualBand_!=0)
			actualLambda=actualBand_;
        const volScalarField& kappa = dort_.kappaLambda(actualLambda);
		const volScalarField& sigma = dort_.sigmaLambda(actualLambda);

        const surfaceScalarField Ji(dAve_ & mesh_.Sf());

        fvScalarMatrix IiEq =
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
          + fvm::Sp((kappa+sigma)*omega_, ILambda_[lambdaI])
          == sigma*inScatter_[lambdaI]
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
		if(actualBand_!=0)
			break;
    }

    return maxResidual;
}

void Foam::radiation::heteroExtDO::setBandScatter
 (const label lambdaI, const volScalarField value)
{
    inScatter_[lambdaI]=value;
}

void Foam::radiation::heteroExtDO::changeBand
 (const label lambdaI)
{
	actualBand_=lambdaI;
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
