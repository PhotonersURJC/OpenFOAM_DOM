/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::radiation::isotropicRadSource

Description

    Isotropic emmision source. Light is equally distributed in all discrete
	ordinates with their central direction leaving the face. Quadrature keeps
	the mesh coordinates system.

SourceFiles
    isotropicRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceisotropicRadSource_H
#define radSourceisotropicRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "booleanSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class isotropicRadSource Declaration
\*---------------------------------------------------------------------------*/

class isotropicRadSource
:
    public radSource
{
	// Private data
		PtrList<booleanSourceStorage> sourceValues_;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        isotropicRadSource(const isotropicRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const isotropicRadSource&);


public:

    //- Runtime type information
    TypeName("isotropicRadSource");


    // Constructors

        //- Construct from components
        isotropicRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        isotropicRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~isotropicRadSource();


    // Member functions

	    //- Return main axis to draw the new quadrature
	    const List<vector> getAxis(label nTheta, label nPhi);
		
            //- Store values of emission
            virtual void generate(const List<vector> directions);

            //- Read radiation properties dictionary
            bool read();

        // Access
		
			//- Acces to specific stored value
            scalar calculate(const word patchID, const label face, const label direction, const label band) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
