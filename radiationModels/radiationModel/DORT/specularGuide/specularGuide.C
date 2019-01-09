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

#include "specularGuide.H"
#include "fvm.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::specularGuide::specularGuide
(
    const Foam::fvPatch& patch,
    const Foam::List<Foam::vector> mainDirs,
	const List<List< label > > facesNeighbours,
	const List<scalar> faceMagnitudes,
	const List<vector> quadVertex,
	const List<face> quadFaces,
	const bool secondOrder,
	const label nLambda
)
:
    patch_(patch),
    indexes_(patch_.size(),(short)(-1)),
	interpolators_(0),
	storedIntensities_(0)
{
	//stores a specularIndexer for each face direction
	const vectorField n(patch_.nf());
	for(label lambdaI = 0; lambdaI < nLambda; lambdaI++)
	{
		PtrList<scalarField> bandIntensities(mainDirs.size());
		forAll(bandIntensities, dirI)
		{
			bandIntensities.set(dirI, new scalarField(patch_.size(),0.0));
		}
		storedIntensities_.append(bandIntensities);
	}
	forAll(patch_, faceI)
	{
		if(indexes_[faceI]==-1)
		{
			vector currentN = n[faceI];
			indexes_[faceI]=interpolators_.size();
			interpolators_.setSize(interpolators_.size()+1);
			interpolators_.set
			(
				interpolators_.size()-1,
				new specularIndexer
				(
					facesNeighbours,
					mainDirs,
					faceMagnitudes,
					currentN,
					quadVertex,
					quadFaces,
					secondOrder
					
				)
			);
			for(label i=faceI+1;i<patch_.size();i++)
			{
				if(n[i]==currentN) //may be introduced a relaxed criteria
				{
					indexes_[i]=interpolators_.size()-1;
				}
			}
		}
	}
}

// * * * * * * * * * * * * * * Member functions * * * * * * * * * * * * * * * //

void Foam::radiation::specularGuide::updateSpecular
(
	const Foam::PtrList<Foam::scalarField> patchIntensity, const label lambdaI
)
{
	
	forAll(patch_, faceId)
	{
		//generates a scalarField for every face
		scalarField faceIntensities(patchIntensity.size());
		forAll(faceIntensities, dirI)
		{
			faceIntensities[dirI] = patchIntensity[dirI][faceId];
		}
		scalarField reflectedIntensities = 
			interpolators_[indexes_[faceId]].getSpec(faceIntensities);
		forAll(reflectedIntensities, dirI)
		{
			storedIntensities_[lambdaI][dirI][faceId] = reflectedIntensities[dirI];
		}
	}
}

Foam::scalarField Foam::radiation::specularGuide::getSpec
(
	const Foam::label rayId, const Foam::label lambdaId
) const
{
	return storedIntensities_[lambdaId][rayId];
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::specularGuide::~specularGuide()
{}

// ************************************************************************* //
