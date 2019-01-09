    
	const wordList solidNames(rp["solid"]);

    PtrList<fvMesh> solidRegions(solidNames.size());
	PtrList<volScalarField> solidTemperatures(solidNames.size());
    forAll(solidNames, i)
    {
        Info<< "Create solid mesh for region " << solidNames[i]
            << " for time = " << runTime.timeName() << nl << endl;
        solidRegions.set
        (
            i,
            new fvMesh
            (
                IOobject
                (
                    solidNames[i],
                    runTime.timeName(),
                    runTime,
                    IOobject::MUST_READ
                )
            )
        );
    }
	    // Initialise solid field pointer lists
    PtrList<radiation::newRadiationModel> solidRadiation(solidRegions.size());
    // Populate solid field pointer lists
    forAll(solidRegions, i)
    {
        Info<< "*** Reading solid radiation properties for region "
            << solidRegions[i].name() << nl << endl;
		// Create a null temperature field to inactivate thermal calculations
		solidTemperatures.set(
		i,
		new volScalarField
		(
		   IOobject
		   (
				"T",
				runTime.timeName(),
				solidRegions[i],
				IOobject::MUST_READ,
				IOobject::AUTO_WRITE
			),
		   solidRegions[i]
		)
		);
        solidRadiation.set
        (
            i,
            radiation::newRadiationModel::New(solidTemperatures[i])
        );
    }

	