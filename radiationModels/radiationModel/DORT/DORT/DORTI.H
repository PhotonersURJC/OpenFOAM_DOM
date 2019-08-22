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

inline Foam::label Foam::radiation::DORT::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::DORT::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::DORT::nRay() const
{
    return 8*nPhi_*nTheta_;
}


inline Foam::label Foam::radiation::DORT::nLambda() const
{
    return nLambda_;
}

inline Foam::label Foam::radiation::DORT::nQuad() const
{
    return nQuad_;
}

inline const Foam::volScalarField& Foam::radiation::DORT::kappa() const
{
    return kappa_;
}


inline const Foam::volScalarField& Foam::radiation::DORT::kappaLambda
(
    const label lambdaI
) const
{
    return kappaLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::sigma() const
{
    return sigma_;
}


inline const Foam::volScalarField& Foam::radiation::DORT::sigmaLambda
(
    const label lambdaI
) const
{
    return sigmaLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::Glambda
(
	const label lambdaI
) const
{
    return Glambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::DORT::QrLambda
(
	const label lambdaI
) const
{
    return QrLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::QinLambda
(
	const label lambdaI
) const
{
    return QinLambda_[lambdaI];
}

inline const Foam::volScalarField& Foam::radiation::DORT::Elambda
(
	const label lambdaI
) const
{
    return Elambda_[lambdaI];
}

inline Foam::label Foam::radiation::DORT::maxIter() const
{
    return maxIter_;
}

inline Foam::scalar Foam::radiation::DORT::convergence() const
{
    return convergence_;
}


// ************************************************************************* //
