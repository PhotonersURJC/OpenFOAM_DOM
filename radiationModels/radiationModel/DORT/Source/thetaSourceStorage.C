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

#include "thetaSourceStorage.H"
#include "constants.H"


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::thetaSourceStorage::thetaSourceStorage(const label nBands, const label nTheta, const label nPhi)
:
    
    Qlist_(nBands),
    nTheta_(nTheta/2),
	nPhi_(nPhi)
{
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //

Foam::radiation::thetaSourceStorage::~thetaSourceStorage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::thetaSourceStorage::setLambertQs(const List<scalar> Qlist, const scalar expFactor, const scalar correctionFactor)
{
	List<scalar> thetaFactors(nTheta_);
	scalar deltaTheta = (Foam::constant::mathematical::pi*0.5)/nTheta_;
	Info << "nTheta: " << nTheta_ << endl;
	Info << "deltaTheta: " << deltaTheta << endl;
	forAll(thetaFactors, thetaI)
	{
		thetaFactors[thetaI] = correctionFactor*Foam::pow(Foam::cos((thetaI+0.5)*deltaTheta),expFactor);
		Info << "thetaFactors[" << thetaI << "]: " << thetaFactors[thetaI] << endl;
	}
	forAll(Qlist_,lambdaI)
	{
		List<scalar> bandQs(nTheta_);
		forAll(thetaFactors, thetaI)
		{
			bandQs[thetaI]=Qlist[lambdaI]*thetaFactors[thetaI];
		}
		Qlist_[lambdaI]=bandQs;
	}
}

Foam::scalar Foam::radiation::thetaSourceStorage::getValue
(
	const label face, const label direction, const label band
) const
{
	label thetaI = direction/(4*nPhi_);
	if(thetaI<nTheta_)
		return Qlist_[band][thetaI];
	return 0.0;
}


// ************************************************************************* //
