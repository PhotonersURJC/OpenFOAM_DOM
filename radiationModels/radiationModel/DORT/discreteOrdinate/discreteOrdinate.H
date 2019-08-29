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
    Foam::radiation::discreteOrdinate

Description
    Abstract class for radiative intensity storage in a given solid angle.
	May be:
		- cleanDiscOrd: Transparent media
		- homoDiscOrd: Absorptive and scattering media, homogeneous coeffs
		- fullDiscOrd: same as homoDiscOrd, but heterogeneous coeffs

SourceFiles
    discreteOrdinate.C

\*---------------------------------------------------------------------------*/

#ifndef discreteOrdinate_H
#define discreteOrdinate_H

#include "blackBodyEmissionRev.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

// Forward declaration of classes
class DORT;
class quadrature;


/*---------------------------------------------------------------------------*\
                    Class discreteOrdinate Declaration
\*---------------------------------------------------------------------------*/

class discreteOrdinate
{
public:

    static const word intensityPrefix;


protected:

    // Private data
	
        //- Refence to the owner DORT object
        const DORT& dort_;

	//- Reference to the owner quadrature
	const quadrature& quad_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Direction
        vector d_;

        //- Average direction vector inside the solid angle
        vector dAve_;

        //- Theta angle
        scalar theta_;

        //- Phi angle
        scalar phi_;

        //- Solid angle
        scalar omega_;

        //- Number of wavelengths/bands
        label nLambda_;

        //- List of pointers to radiative intensity fields for given wavelengths
        PtrList<volScalarField> ILambda_;
		
        //- Global ray id - incremented in constructor
        static label rayId;

        //- My ray Id
        label myRayId_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        discreteOrdinate(const discreteOrdinate&);

        //- Disallow default bitwise assignment
        void operator=(const discreteOrdinate&);


public:

    // Constructors

        //- Construct form components
        discreteOrdinate
        (
            const DORT& dort,
	    const quadrature& quad,
            const fvMesh& mesh,
            const scalar phi,
            const scalar theta,
            const scalar deltaPhi,
            const scalar deltaTheta,
            const label lambda,
            const label rayId,
	    const List<vector> mainAxis
        );


    //- Destructor
    virtual ~discreteOrdinate();


    // Member functions

        // Edit

        //- Update radiative intensity on i direction
        virtual scalar correct() = 0;

	    virtual void setBandScatter (const label lambdaI, const volScalarField value) = 0;
		
		virtual void changeBand (const label lambdaI) = 0;

        // Access

            //- Return direction
            inline const vector& d() const;

            //- Return the average vector inside the solid angle
            inline const vector& dAve() const;

            //- Return the number of bands
            inline scalar nLambda() const;

            //- Return the phi angle
            inline scalar phi() const;

            //- Return the theta angle
            inline scalar theta() const;

            //- Return the solid angle
            inline scalar omega() const;

            //- Return the radiative intensity for a given wavelength
            inline const volScalarField& ILambda(const label lambdaI) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace radiation
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "discreteOrdinateI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
