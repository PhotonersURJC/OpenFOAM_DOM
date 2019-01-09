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

#include "specularMatrix.H"
#include "fvm.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::specularMatrix::specularMatrix
(
	const PrimitivePatch<face,List,pointField,point> quadraturePatch,
	const PrimitivePatch<face,List,pointField,point> reflectedPatch,
	const List<vector> mainDirs,
	const vector faceN
)
:
fromPatch_(reflectedPatch),
toPatch_(quadraturePatch),
interpolator_(fromPatch_,toPatch_),
dirOut_(mainDirs.size(),false)
{
	forAll(mainDirs, dirI)
	{
		if((-faceN & mainDirs[dirI]) > 0.0)
		{
			dirOut_[dirI]=true;
		}
	}
}

Foam::scalarField Foam::radiation::specularMatrix::getSpec
(
	const Foam::scalarField originalI
) const
{
	return interpolator_.faceInterpolate(originalI);
}

bool Foam::radiation::specularMatrix::dirOut
(
	const label dirI
) const
{
	return dirOut_[dirI];
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::specularMatrix::~specularMatrix()
{}

// ************************************************************************* //
