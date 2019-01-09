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

#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::booleanSourceStorage::booleanSourceStorage(const label nBands, const label nFaces, const label nDirections)
:
    
    Qlist_(nBands),
    criteria_(nFaces)
{
	forAll(criteria_, critI)
	{
		List<bool> currentCrit(nDirections,false);
		criteria_[critI]=currentCrit;
	}
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::booleanSourceStorage::~booleanSourceStorage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::booleanSourceStorage::setQs(const List<scalar> Qlist)
{
	forAll(Qlist_,lambdaI)
	{
		Qlist_[lambdaI]=Qlist[lambdaI];
	}
}

void Foam::radiation::booleanSourceStorage::setBools
(
	const label face, const label direction, const bool criteria
)
{
	if(criteria)
	{
		criteria_[face][direction]=true;
	}
}

Foam::scalar Foam::radiation::booleanSourceStorage::getValue
(
	const label face, const label direction, const label band
) const
{
	if(criteria_[face][direction])
	{
		return Qlist_[band];
	}
	return 0.0;
}


// ************************************************************************* //
