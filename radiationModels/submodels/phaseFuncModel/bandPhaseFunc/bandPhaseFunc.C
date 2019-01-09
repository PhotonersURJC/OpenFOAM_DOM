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

#include "bandPhaseFunc.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(bandPhaseFunc, 0);

        addToRunTimeSelectionTable
        (
            phaseFuncModel,
            bandPhaseFunc,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::bandPhaseFunc::bandPhaseFunc
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    phaseFuncModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Model")),
    phaseFunctions_(0)
{
	wordList bandNames = coeffsDict_.toc();
	phaseFunctions_.resize(bandNames.size());
	forAll(bandNames, lambdaI)
	{
		dictionary currentBandDict = coeffsDict_.subDict(bandNames[lambdaI]);
		dictionary newBandDict = new dictionary();
		word modelName = currentBandDict.lookupOrDefault<word>("model", "isotropicPhaseFunc");
		newBandDict.add(keyType("phaseFunction"),modelName);
		modelName += "Coeffs";
		newBandDict.add(keyType(modelName),currentBandDict.subDict("modelCoeffs"));
		phaseFunctions_.set
		(
			lambdaI,
			phaseFuncModel::New(newBandDict, mesh_).ptr()
		);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::bandPhaseFunc::~bandPhaseFunc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::bandPhaseFunc::MatrixGen
(const List<vector> directions, const label model)
{
	forAll(phaseFunctions_, lambdaI)
	{
		phaseFunctions_[lambdaI].MatrixGen(directions,model);
	}
}

Foam::scalar Foam::radiation::bandPhaseFunc::sctrFactor
(const label model, const label band, const label directionFrom, const label directionTo)
const
{
	return phaseFunctions_[band].sctrFactor(model, band, directionFrom, directionTo);
}


// ************************************************************************* //
