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
    Struct to storage source emmsion power in all directions in a boundary.
	It stores a single double value in each face, and a boolean value
	in each face and direction, identifying directions where there is emmision

Class
    Foam::radiation::booleanSourceStorage

Description
    

SourceFile
    booleanSourceStorage.C

\*---------------------------------------------------------------------------*/

#ifndef radiationModelBooleanSourceStorage_H
#define radiationModelBooleanSourceStorage_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "IOdictionary.H"

namespace Foam
{

namespace radiation
{

/*---------------------------------------------------------------------------*\
                    Class booleanSourceStorage Declaration
\*---------------------------------------------------------------------------*/

class booleanSourceStorage
{
protected:

    // Protected data

	//- patch Q in each Band
	List<scalar> Qlist_;
	
	//- Source emmision in each face and direction
	List<List<bool> > criteria_;

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        booleanSourceStorage(const booleanSourceStorage&);

        //- Disallow default bitwise assignment
        void operator=(const booleanSourceStorage&);


public:

    // Constructors

        //- Construct from components
        booleanSourceStorage(const label nBands, const label nFaces, const label nDirections);


    //- Destructor
    ~booleanSourceStorage();


    // Member Functions

        // Edit

            //- Store Fluxes
            void setQs(const List<scalar> Qlist);


			//- Store boolean values
			void setBools(const label face, const label direction, const bool criteria);


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
