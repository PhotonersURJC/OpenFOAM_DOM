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

#ifndef radiationSpecularMatrix_H
#define radiationSpecularMatrix_H

#include "volFields.H"
#include "PatchToPatchInterpolation.H"
//#include "MeshedSurface.H"
#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularIndexer Declaration
\*---------------------------------------------------------------------------*/

class specularMatrix
{
				
		PrimitivePatch<face,List,vectorField,point> fromPatch_;
		
		PrimitivePatch<face,List,vectorField,point> toPatch_;
		
		PatchToPatchInterpolation<PrimitivePatch<face,List,vectorField,point>,PrimitivePatch<face,List,vectorField,point> > interpolator_;
		
		List<bool> dirOut_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularMatrix(const specularMatrix&);

        //- Disallow default bitwise assignment
        void operator=(const specularMatrix&);


public:

    // Constructors
		
        //- Construct from components
        specularMatrix
        (
			const PrimitivePatch<face,List,pointField,point> quadraturePatch,
			const PrimitivePatch<face,List,pointField,point> reflectedPatch,
			const List<vector> mainDirs,
			const vector faceN
        );

    //- Destructor
    ~specularMatrix();


    // Member functions

	    //- Return the specular index for a given direction
            scalarField getSpec(const scalarField originalI) const;
			
		//- Determines if a direction is leaving the face, or not
			bool dirOut(const label dirI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
