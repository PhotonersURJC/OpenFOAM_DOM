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
    Foam::radiation::specularGuide

Description
    Class to storage specular indexes for a given face in a region.
	Stores a specularIndexer per face direction (a flat face stores only
	a single specularIndexer).

SourceFiles
    specularGuide.C

\*---------------------------------------------------------------------------*/

#ifndef radiationSpecularGuide_H
#define radiationSpecularGuide_H

#include "specularIndexer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{


/*---------------------------------------------------------------------------*\
                    Class specularGuide Declaration
\*---------------------------------------------------------------------------*/

class specularGuide
{

        //- Reference to the patch
        const fvPatch& patch_;
		
		List<short> indexes_;

		PtrList<specularIndexer> interpolators_;
		
		//- Intensity in [band][direction][face]
        List<PtrList<scalarField> > storedIntensities_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        specularGuide(const specularGuide&);

        //- Disallow default bitwise assignment
        void operator=(const specularGuide&);

public:

    // Constructors

        //- Construct from components
        specularGuide
        (
            const Foam::fvPatch& patch,
			const Foam::List<Foam::vector> mainDirs,
			const List<List< label > > facesNeighbours,
			const List<scalar> faceMagnitudes,
			const List<vector> quadVertex,
			const List<face> quadFaces,
			const bool secondOrder,
			const label nLambda
        );


    //- Destructor
    ~specularGuide();


    // Member functions
		//- Assignment
			void updateSpecular(const PtrList<scalarField> patchIntensity, const label lambdaI);

	    //- Return the reflecteds Intensities for a given face
            scalarField getSpec(const label rayId, const label lambdaId) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
