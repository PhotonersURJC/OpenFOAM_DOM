/*---------------------------------------------------------------------------*\
             Discrete Ordinate Method Radiation Model for OpenFOAM.


Code corresponding to the article entitled:


"Improved Discrete Ordinate Method for accurate simulation radiation transport
                  using solar and LED light sources"


by:


José Moreno, Cintia Casado, Javier Marugán

Department of Chemical and Environmental Technology,

ESCET, Universidad Rey Juan Carlos,

C/Tulipán s/n, 28933 Móstoles (Madrid), Spain

Tel. +34 91 664 7466; E-mail: javier.marugan@urjc.es

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
