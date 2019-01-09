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

#include "discreteOrdinate.H"
#include "fvm.H"
#include "DORT.H"
#include "quadrature.H"
#include "constants.H"

using namespace Foam::constant;


const Foam::word
Foam::radiation::discreteOrdinate::intensityPrefix("ILambda");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::discreteOrdinate::discreteOrdinate
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
    dort_(dort),
    quad_(quad),
    mesh_(mesh),
    d_(vector::zero),
    dAve_(vector::zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    myRayId_(rayId)
{

    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    if(omega_<0) { omega_ = -omega_;}
    d_ = vector(sinTheta*sinPhi*mainAxis[0])+ (sinTheta*cosPhi*mainAxis[1]) + (cosTheta*mainAxis[2]);
    dAve_ = omega_ * d_;
    autoPtr<volScalarField> IDefaultPtr;
    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(quad.modelNumber()) + "_" + name(rayId) + "_" +  name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );
		
        // check if field exists and can be read
        if (IHeader.typeHeaderOk<volScalarField>(true))
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IBand field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IBand_" + name(lambdaI),
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }
			
            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::discreteOrdinate::~discreteOrdinate()
{}



// ************************************************************************* //
