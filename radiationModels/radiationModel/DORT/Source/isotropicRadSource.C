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

#include "isotropicRadSource.H"
#include "phaseFuncModel.H"
#include "radExtintionModel.H"
#include "constants.H"
#include "fvm.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(isotropicRadSource, 0);
        addToSourceRunTimeSelectionTables(isotropicRadSource);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::radiation::isotropicRadSource::initialise()
{
	sourceValues_.setSize(patchesID_.size());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::isotropicRadSource::isotropicRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands)
:
    radSource(typeName, mesh, sourceIndex, nBands),
	sourceValues_(0)

{
    initialise();
}


Foam::radiation::isotropicRadSource::isotropicRadSource
(
    const dictionary& dict,
    const fvMesh& mesh,
    const label sourceIndex,
	const label nBands
)
:
    radSource(typeName, dict, mesh, sourceIndex, nBands),
	sourceValues_(0)
{
    initialise();
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::isotropicRadSource::~isotropicRadSource()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::radiation::isotropicRadSource::read()
{
    if (radSource::read())
    {


        return true;
    }
    else
    {
        return false;
    }
}

const Foam::List<Foam::vector> Foam::radiation::isotropicRadSource::getAxis
(Foam::label nTheta, Foam::label nPhi)
{
    List<vector> mainAxis(3,vector());
    vector Vx(1,0,0); vector Vy(0,1,0); vector Vz(0,0,1);
    mainAxis[0] = Vx; mainAxis[1] = Vy; mainAxis[2] = Vz;
    return mainAxis;
}

Foam::scalar Foam::radiation::isotropicRadSource::calculate
(
	const word patchID, const label face, const label direction, const label band
) const
{
	forAll(patchesID_, patch)
	{
		if(patchID == patchesID_[patch])
		{
			return sourceValues_[patch].getValue(face, direction, band);
		}
	}
	return 0.0;
}

void Foam::radiation::isotropicRadSource::generate(const List<vector> directions)
{
    label patchIndex = 0;
    while (patchIndex < patchesID_.size())
    {
        forAll(mesh_.boundary(), patchI)
        {
	    const polyPatch currentPatch(mesh_.boundary()[patchI].patch());
            if(currentPatch.name()==patchesID_[patchIndex])
            {
				sourceValues_.set(patchIndex,new booleanSourceStorage(nBands_,currentPatch.size(),directions.size()));
				List<scalar> modQs = patchesQ_[patchIndex];
				forAll(modQs, currentQ)
				{
					modQs[currentQ]=modQs[currentQ]/constant::mathematical::pi;
				}
				sourceValues_[patchIndex].setQs(modQs);
				forAll(currentPatch, faceI)
				{
					forAll(directions, currentDir)
					{
						sourceValues_[patchIndex].setBools(faceI,currentDir,true);
					}
				}
			}
        }
        patchIndex++;
    }
}



// ************************************************************************* //
