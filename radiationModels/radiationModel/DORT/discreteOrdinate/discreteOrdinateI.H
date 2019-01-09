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

inline const Foam::vector& Foam::radiation::discreteOrdinate::d() const
{
    return d_;
}


inline const Foam::vector& Foam::radiation::discreteOrdinate::dAve() const
{
    return dAve_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::nLambda() const
{
    return nLambda_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::phi() const
{
    return phi_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::theta() const
{
    return theta_;
}


inline Foam::scalar Foam::radiation::discreteOrdinate::omega() const
{
    return omega_;
}

inline const Foam::volScalarField&
Foam::radiation::discreteOrdinate::ILambda
(
    const label lambdaI
) const
{
    return ILambda_[lambdaI];
}


// ************************************************************************* //
