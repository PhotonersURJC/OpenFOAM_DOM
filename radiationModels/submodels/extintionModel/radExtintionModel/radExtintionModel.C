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

#include "radExtintionModel.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(radExtintionModel, 0);					
        defineRunTimeSelectionTable(radExtintionModel, dictionary);		
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radExtintionModel::radExtintionModel
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    dict_(dict),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::radExtintionModel::~radExtintionModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::radExtintionModel::kappa(const label bandI) const
{
    return kappaDisp(bandI) + kappaCont(bandI);
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::radExtintionModel::sigma(const label bandI) const
{
    return sigmaDisp(bandI) + sigmaCont(bandI);
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::radExtintionModel::kappaCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "kappaCont",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::radExtintionModel::kappaDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "kappaDisp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}

Foam::volScalarField Foam::radiation::radExtintionModel::fullKappa() const
{
	volScalarField result = this->kappa(0);
	for(int lambdaI = 1; lambdaI < this->nBands(); lambdaI++)
		result+= this->kappa(lambdaI);
	return result;
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::radExtintionModel::sigmaCont(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaCont",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::radExtintionModel::sigmaDisp(const label bandI) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "sigmaDisp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar("zero", dimless/dimLength, 0.0)
        )
    );
}


Foam::label Foam::radiation::radExtintionModel::nBands() const
{
    return pTraits<label>::one;
}


bool Foam::radiation::radExtintionModel::isGrey() const
{
    return false;
}

bool Foam::radiation::radExtintionModel::isConstant() const
{
    return false;
}

bool Foam::radiation::radExtintionModel::isHomogeneous() const
{
    return false;
}

void Foam::radiation::radExtintionModel::kappaCorrect
(
    volScalarField& kappa,
    PtrList<volScalarField>& kappaj
) const
{
    kappa = this->kappa();
    kappaj[0] =  kappa;
}

void Foam::radiation::radExtintionModel::sigmaCorrect
(
    volScalarField& sigma,
    PtrList<volScalarField>& sigmaj
) const
{
    sigma = this->sigma();
    sigmaj[0] =  sigma;
}

// ************************************************************************* //
