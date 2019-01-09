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

#include "newRadiationModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


Foam::autoPtr<Foam::radiation::newRadiationModel>
Foam::radiation::newRadiationModel::New
(
    const volScalarField& T
)
{
    IOobject radIO
    (
        "radiationProperties",
        T.time().constant(),
        T.mesh(),
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE,
        false
    );
    word modelType("none");
    if (radIO.good())
    {
        IOdictionary(radIO).lookup("radiationModel") >> modelType;
    }
    else
    {
        Info<< "Radiation model not active: radiationProperties not found"
            << endl;
    }
    Info<< "Selecting radiationModel " << modelType << endl;
    TConstructorTable::iterator cstrIter =
        TConstructorTablePtr_->find(modelType);

    if (cstrIter == TConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "noTradiationModel::New(const fvMesh&)"
        )   << "Unknown noTradiationModel type "
            << modelType << nl << nl
            << "Valid noTradiationModel types are:" << nl
            << TConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }
    return autoPtr<newRadiationModel>(cstrIter()(T));
}

Foam::autoPtr<Foam::radiation::newRadiationModel>
Foam::radiation::newRadiationModel::New
(
    const dictionary& dict,
    const volScalarField& T
)
{
    const word modelType(dict.lookup("radiationModel"));

    Info<< "Selecting radiationModel " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "noTradiationModel::New(const dictionary&, const fvMesh&)"
        )   << "Unknown radiationModel type "
            << modelType << nl << nl
            << "Valid radiationModel types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<newRadiationModel>(cstrIter()(dict, T));
}


// ************************************************************************* //
