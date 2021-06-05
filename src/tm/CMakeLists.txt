if (USE_CEMHYD)
    set (tm_cemhyd
        Materials/cemhyd/cemhydmat.C
        )
endif ()

set (tm_emodel
    EngineeringModels/nltransienttransportproblem.C
    EngineeringModels/nonstationarytransportproblem.C
    EngineeringModels/stationarytransportproblem.C
    EngineeringModels/transienttransportproblem.C
    )

set (tm_element
<<<<<<< HEAD
    transportelement.C
    quad1_ht.C
    tr1_ht.C
    quadaxisym1_ht.C
    traxisym1_ht.C
    brick1_ht.C
    qbrick1_ht.C
    wedge_ht.C
    qwedge_ht.C
    tetrah1_ht.C
    tr1darcy.C
    latticetransportelement.C
    lattice2d_mt.C
    )

set (tm_boundary
    gravitypressure.C
=======
    Elements/transportelement.C
    Elements/quad1_ht.C
    Elements/qquad1_ht.C
    Elements/tr1_ht.C
    Elements/quadaxisym1_ht.C
    Elements/traxisym1_ht.C
    Elements/brick1_ht.C
    Elements/qbrick1_ht.C
    Elements/wedge_ht.C
    Elements/qwedge_ht.C
    Elements/tetrah1_ht.C
    Elements/tr1darcy.C
    Elements/LatticeElements/latticetransportelement.C
    Elements/LatticeElements/lattice2d_mt.C
    Elements/LatticeElements/lattice3d_mt.C
    )

set (tm_boundary
    BoundaryCondition/gravitypressure.C
    BoundaryCondition/transportgradientneumann.C
    BoundaryCondition/transportgradientdirichlet.C
    BoundaryCondition/transportgradientperiodic.C
>>>>>>> bp2/master
    )

set (tm_cs
    simpletransportcrosssection.C
    )

set (tm_material
<<<<<<< HEAD
    hydram.C
    hydratingisoheatmat.C
    hydratinghemomat.C
    hydratingconcretemat.C
    hemotkmat.C
    hemokunzelmat.C
    hemobaznajmat.C
    isoheatmat.C
    isoheatadvmat.C
    isomoisturemat.C
    isolinmoisturemat.C
    bazantnajjarmat.C
    nlisomoisturemat.C
    nonlinearmassmat.C
    anisomassmat.C
    rvestokesflow.C
=======
    Materials/anisomassmat.C
    Materials/bazantnajjarmat.C
    Materials/hemobaznajmat.C
    Materials/hemokunzelmat.C
    Materials/hemotkmat.C
    Materials/hydram.C
    Materials/hydratingconcretemat.C
    Materials/hydratinghemomat.C
    Materials/hydratingisoheatmat.C
    Materials/isoheatmat.C
    Materials/isolinmoisturemat.C
    Materials/isomoisturemat.C
    Materials/LatticeMaterials/latticetransmat.C
    Materials/nlisomoisturemat.C
    Materials/nonlinearmassmat.C
    Materials/rvestokesflow.C
    Materials/transportmaterial.C
    Materials/twophasemat.C
>>>>>>> bp2/master
    )

set (tm
    ${tm_emodel}
    ${tm_element}
    ${tm_material}
    ${tm_cs}
    ${tm_new}
    ${tm_cemhyd}
    ${tm_hydram}
    ${tm_boundary}
    )

add_library (tm
    ${LIB_TYPE}
    ${tm}
    )
