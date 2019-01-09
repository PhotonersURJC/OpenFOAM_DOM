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

#include "coeffPhaseFunc.H"
#include "addToRunTimeSelectionTable.H"
#include "constants.H"

using namespace Foam::constant;
using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(coeffPhaseFunc, 0);

        addToRunTimeSelectionTable
        (
            phaseFuncModel,
            coeffPhaseFunc,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::coeffPhaseFunc::coeffPhaseFunc
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    phaseFuncModel(dict, mesh),
    coeffsDict_(dict.subDict(typeName + "Model")),
    phaseFunctions_(0)
{
	word modelName = coeffsDict_.lookupOrDefault<word>("model", "isotropicPhaseFunc");
	label nBands = coeffsDict_.lookupOrDefault<label>("nBands", 1);
	wordList coeffsNames = coeffsDict_.toc();
	phaseFunctions_.resize(nBands);
	forAll(phaseFunctions_,lambdaI)
	{
		dictionary newBandDict = new dictionary();
		newBandDict.add(keyType("phaseFunction"),modelName);
		word modelCoeffs = modelName + "Coeffs";
		if(coeffsNames.size()>2)
		{
			dictionary currentBandCoeffs = new dictionary();
			for(label coeff = 2; coeff < coeffsNames.size(); coeff++)
			{
				word currentBandName = "band"+name(lambdaI);
				currentBandCoeffs.add(coeffsDict_.subDict(coeffsNames[coeff]).lookupEntry(currentBandName,false,true));
				currentBandCoeffs.changeKeyword(currentBandName,coeffsNames[coeff]);
			}
			newBandDict.add(keyType(modelCoeffs),currentBandCoeffs);
		}
		phaseFunctions_.set
		(
			lambdaI,
			phaseFuncModel::New(newBandDict, mesh_).ptr()
		);
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::coeffPhaseFunc::~coeffPhaseFunc()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::coeffPhaseFunc::MatrixGen
(const List<vector> directions, const label model)
{
	forAll(phaseFunctions_, lambdaI)
	{
		phaseFunctions_[lambdaI].MatrixGen(directions,model);
	}
}

Foam::scalar Foam::radiation::coeffPhaseFunc::sctrFactor
(const label model, const label band, const label directionFrom, const label directionTo)
const
{
	return phaseFunctions_[band].sctrFactor(model, band, directionFrom, directionTo);
}


// ************************************************************************* //
