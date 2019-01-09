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
