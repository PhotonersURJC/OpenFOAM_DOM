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
