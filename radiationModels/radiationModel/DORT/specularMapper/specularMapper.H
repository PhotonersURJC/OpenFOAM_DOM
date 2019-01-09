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
    Foam::radiation::discreteOrdinate

Description
    Radiation intensity for a ray in a given direction

SourceFiles
    discreteOrdinate.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularIndexer_H
#define radiationSpecularIndexer_H

//#include "fvPatch.H"
//#include "fvMatrices.H"
#include "specularMatrix.H"

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

        //- Reference to the patch
        const fvPatch& patch_;
		
		List<short> indexes_;
		
        //- Intensity in [band][direction][face]
        List<PtrList<scalarField> > storedIntensities_;

		PtrList<specularMatrix> interpolators_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularIndexer(const specularIndexer& /*toCopy*/);

        //- Disallow default bitwise assignment
        void operator=(const specularIndexer&);

public:

    // Constructors
	
		//- Null constructor
		/*specularIndexer();*/

        //- Construct from components
        specularIndexer
        (
            const fvPatch& patch,
            const List<vector> mainDirs,
			const List<vector> vertex,
			const List<face> faces,
			const label nLambda
        );


    //- Destructor
    ~specularIndexer();


    // Member functions
		//- Assignment
			void updateSpecular(const PtrList<scalarField> patchIntensity, const label lambdaI);

	    //- Return the specular index for a given face
            scalarField getSpec(const label rayId, const label lambdaId) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
