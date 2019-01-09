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

#include "sctrMatrix.H"
#include "fvm.H"
#include "constants.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::sctrMatrix::sctrMatrix
(
    const Foam::label directions
)
:
    sctrFactors_(directions*directions),
	directions_(directions)
{

}

Foam::scalar& Foam::radiation::sctrMatrix::Item
(const label i, const label j)
{
	return sctrFactors_[i*directions_+j];
}

Foam::scalar Foam::radiation::sctrMatrix::constItem
(const label i, const label j) const
{
	return sctrFactors_[i*directions_+j];
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::sctrMatrix::~sctrMatrix()
{}

// ************************************************************************* //
