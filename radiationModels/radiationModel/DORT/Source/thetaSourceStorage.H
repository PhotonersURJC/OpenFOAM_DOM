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

Namespace
    Foam::radiation

Description
    Struct to storage source emmision power in all directions in a boundary.
	It stores a double value in each band and theta value of the quadrature.
	Every face in the boundary emits exactly the same intensities
Class
    Foam::radiation::thetaSourceStorage

Description
    

SourceFile
    thetaSourceStorage.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelThetaSourceStorage_H
#define radiationModelThetaSourceStorage_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOdictionary.H"

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class thetaSourceStorage Declaration
\*---------------------------------------------------------------------------*/

class thetaSourceStorage
{
protected:

    // Protected data

	//- patch Q in each Band and Theta Index
	List<List<scalar> > Qlist_;
	
	//quadrature theta divisions
	label nTheta_;
	
	//quadrature phi divisiones
	label nPhi_;
	

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        thetaSourceStorage(const thetaSourceStorage&);

        //- Disallow default bitwise assignment
        void operator=(const thetaSourceStorage&);


public:

    // Constructors

        //- Construct from components
        thetaSourceStorage(const label nBands, const label nTheta, const label nPhi);


    //- Destructor
    ~thetaSourceStorage();


    // Member Functions

        // Edit

            //- Store Fluxes
            void setLambertQs(const List<scalar> Qlist, const scalar expFactor, const scalar correctionFactor);

        // Access

            //- Access to specific value
            scalar getValue(const label face, const label direction, const label band) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
