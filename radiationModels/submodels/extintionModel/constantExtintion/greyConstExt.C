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

#include "greyConstExt.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyConstExt, 0);

        addToRunTimeSelectionTable
        (
            radExtintionModel,
            greyConstExt,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyConstExt::greyConstExt
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    radExtintionModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Coeffs")),
    kappa_(coeffsDict_.lookup("absorptivity")),
	sigma_(coeffsDict_.lookup("scatter"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyConstExt::~greyConstExt()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyConstExt::kappaCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "kappa",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            kappa_
        )
    );

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyConstExt::sigmaCont(const label bandI) const
{
    tmp<volScalarField> ts
    (
        new volScalarField
        (
            IOobject
            (
                "sigma",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            sigma_
        )
    );

    return ts;
}


// ************************************************************************* //
