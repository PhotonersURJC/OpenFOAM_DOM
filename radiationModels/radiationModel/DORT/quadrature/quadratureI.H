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

inline Foam::label 
Foam::radiation::quadrature::modelNumber() const
{
    return modelNumber_;
}

inline const Foam::radiation::discreteOrdinate&
Foam::radiation::quadrature::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::quadrature::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::quadrature::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::quadrature::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::quadrature::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::quadrature::nLambda() const
{
    return nLambda_;
}


inline const Foam::volScalarField& Foam::radiation::quadrature::Glambda
(
	const label lambdaI
) const
{
    return Glambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::quadrature::QrLambda
(
	const label lambdaI
) const
{
    return QrLambda_[lambdaI];
}

inline Foam::volScalarField& Foam::radiation::quadrature::futQrLambda
(
	const label lambdaI
)
{
    return futQrLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::quadrature::QinLambda
(
	const label lambdaI
) const
{
    return QinLambda_[lambdaI];
}

inline Foam::volScalarField& Foam::radiation::quadrature::futQinLambda
(
	const label lambdaI
)
{
    return futQinLambda_[lambdaI];
}

inline Foam::scalar Foam::radiation::quadrature::omegaMax() const
{
    return omegaMax_;
}

inline Foam::scalar Foam::radiation::quadrature::sourceFlux
(const word patchID, const label face, const label rayId, const label band) const
{
    return radSource_.calculate(patchID,face,rayId,band);
}

inline const Foam::List<Foam::vector> Foam::radiation::quadrature::directions() const
{
	return directions_;
}

// ************************************************************************* //
