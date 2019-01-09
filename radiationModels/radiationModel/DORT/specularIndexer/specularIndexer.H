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
Class
    Foam::radiation::specularIndexer

Description
    Storage of specular reflection directions for each discrete Ordinate in
	a given face.

SourceFiles
    specularIndexer.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularIndexer_H
#define radiationSpecularIndexer_H

#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularIndexer Declaration
\*---------------------------------------------------------------------------*/

class specularIndexer
{
		//- Determines where to look for spotters, if any (else -1)
		List<short> outers_;
		
		//- Reflected incoming DOs looking at current outgoing DO
		List<List<short > > spotters_;
		
		//- Factors to add intensities of spotters
		List<List<scalar> > viewFactors_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const specularIndexer&);

        //- Disallow default bitwise assignment
        void operator=(const specularIndexer&);


public:

    // Constructors
		
        //- Construct from components
        specularIndexer
        (
			const List<List< label > > facesNeighbours,
			const List<vector> mainDirs,
			const List<scalar> faceMagnitudes,
			const vector faceN,
			const List<vector> quadVertex,
			const List<face> quadFaces,
			const bool secondOrder
        );

    //- Destructor
    ~specularIndexer();


    // Member functions

	    //- Return the reflected Field for a given face
            scalarField getSpec(const scalarField originalI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
