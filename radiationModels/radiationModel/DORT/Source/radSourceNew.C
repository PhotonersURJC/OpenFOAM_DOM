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

#include "radSource.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::radiation::radSource>
Foam::radiation::radSource::New
(
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
{
    IOobject radIO
    (
        "radSource" + Foam::name(sourceIndex),
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    word modelType("none");
    if (radIO.good())
    {
        IOdictionary(radIO).lookup("sourceType") >> modelType;
    }

    meshConstructorTable::iterator cstrIter =
        meshConstructorTablePtr_->find(modelType);

    if (cstrIter == meshConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "radSource::New(const fvMesh&)"
        )   << "Unknown radSource type "
            << modelType << nl << nl
            << "Valid radSource types are:" << nl
            << meshConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<radSource>(cstrIter()(mesh, sourceIndex, nBands));
}


Foam::autoPtr<Foam::radiation::radSource>
Foam::radiation::radSource::New
(
    const dictionary& dict,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
{
    const word modelType(dict.lookup("sourceType"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "radSource::New(const dictionary&, const fvMesh&)"
        )   << "Unknown radiSource type "
            << modelType << nl << nl
            << "Valid radSource types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<radSource>(cstrIter()(dict, mesh, sourceIndex, nBands));
}


// ************************************************************************* //
