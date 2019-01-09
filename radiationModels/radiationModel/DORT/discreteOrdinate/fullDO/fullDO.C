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

#include "fullDO.H"
#include "fvm.H"
#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::fullDO::fullDO
(
    const DORT& dort,
    const quadrature& quad,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const label rayId,
    const List<vector> mainAxis
)
:
    discreteOrdinate
    (
        dort,
        quad,
        mesh,
        phi,
        theta,
        deltaPhi,
        deltaTheta,
        nLambda,
        rayId,
        mainAxis
    ),
    inScatter_(nLambda),
	actualBand_(0)
{
    forAll(ILambda_, lambdaI)
    {
	inScatter_.set
	(
	    lambdaI,
	    new volScalarField
	    (
		IOobject
		(
		    "inScatter",
		    mesh_.time().timeName(),
		    mesh_,
		    IOobject::NO_READ,
		    IOobject::NO_WRITE
		),
		mesh_,
		dimensionedScalar("inScatter", dimMass/pow3(dimTime), 0.0)
	    )
	);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::fullDO::~fullDO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiation::fullDO::correct
()
{
    scalar maxResidual = -GREAT;

    forAll(ILambda_, lambdaI)
    {
		label actualLambda=lambdaI;
		if(actualBand_!=0)
			actualLambda=actualBand_;
        const volScalarField& kappa = dort_.kappaLambda(actualLambda);
		const volScalarField& sigma = dort_.sigmaLambda(actualLambda);

        const surfaceScalarField Ji(dAve_ & mesh_.Sf());

        fvScalarMatrix IiEq =
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
          + fvm::Sp((kappa+sigma)*omega_, ILambda_[lambdaI])
          == sigma*inScatter_[lambdaI] + dort_.Elambda(actualLambda)
        );
        

        IiEq.relax();
        const solverPerformance ILambdaSol = solve
        (
            IiEq,
            mesh_.solver("Ii")
        );

        const scalar initialRes =
            ILambdaSol.initialResidual()*omega_/quad_.omegaMax();
		ILambda_[lambdaI].max(0.0);
        maxResidual = max(initialRes, maxResidual);
		if(actualBand_!=0)
			break;
    }

    return maxResidual;
}

void Foam::radiation::fullDO::setBandScatter
 (const label lambdaI, const volScalarField value)
{
    inScatter_[lambdaI]=value;
}

void Foam::radiation::fullDO::changeBand
 (const label lambdaI)
{
	actualBand_=lambdaI;
	IOobject IHeader
    (
        "IBand_" +  name(lambdaI),
        mesh_.time().timeName(),
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    );
	if (IHeader.good())
	{
		ILambda_.set
		(
			0,
			new volScalarField(IHeader, mesh_)
		);
	}
}


// ************************************************************************* //
