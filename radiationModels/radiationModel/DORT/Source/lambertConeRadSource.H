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
    Foam::radiation::lambertConeRadSource

Description

    lambertian radiation distribution in an hemisphere. Quadrature rotates
	to fit its Z-axis in the source central direction.

SourceFiles
    lambertConeRadSource.C

\*---------------------------------------------------------------------------*/

#ifndef radSourceLambertConeRadSource_H
#define radSourceLambertConeRadSource_H

#include "radSource.H"
#include "fvMatrices.H"
#include "thetaSourceStorage.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

/*---------------------------------------------------------------------------*\
                           Class lambertConeRadSource Declaration
\*---------------------------------------------------------------------------*/

class lambertConeRadSource
:
    public radSource
{

	PtrList<thetaSourceStorage> sourceValues_;
	scalar expFactor_;
	scalar correctionFactor_;
	

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        lambertConeRadSource(const lambertConeRadSource&);

        //- Disallow default bitwise assignment
        void operator=(const lambertConeRadSource&);


public:

    //- Runtime type information
    TypeName("lambertConeRadSource");


    // Constructors

        //- Construct from components
        lambertConeRadSource(const fvMesh& mesh, const label sourceIndex, const label nBands);

        //- Construct from components
        lambertConeRadSource(const dictionary& dict, const fvMesh& mesh, const label sourceIndex, const label nBands);


    //- Destructor
    virtual ~lambertConeRadSource();


    // Member functions

        // Edit

	    //- Return main axis to draw the new quadrature
	    const List<vector> getAxis(label nTheta, label nPhi);

            //- Store values of emission
            virtual void generate(const List<vector> directions);
			
            //- Read radiation properties dictionary
            bool read();

        // Access
		
			//- Acces to specific stored value
            virtual scalar calculate(const word patchID, const label face, const label direction, const label band) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
