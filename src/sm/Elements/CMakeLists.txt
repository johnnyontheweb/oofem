set (sm_element
    structuralelementevaluator.C
    igaelements.C
    structuralelement.C
    nlstructuralelement.C
    structural2delement.C
    structural3delement.C
    3D/space3delementevaluator.C
    3D/lspace.C 
    3D/lspacebb.C
    3D/qspace.C 
    3D/q27space.C
    3D/ltrspace.C 
    3D/ltrspaceboundary.C
    3D/ltrspaceboundarytruss.C
    3D/ltrspaceboundarymembrane.C
    3D/ltrspaceboundarybeam.C
    3D/ltrspaceboundaryplate.C
    3D/ltrspaceboundaryvoigt.C
    3D/qtrspace.C 
    3D/lwedge.C 
    3D/qwedge.C 
    3D/macrolspace.C
    Bars/truss1d.C 
    Bars/truss2d.C 
    Bars/truss3d.C
    Bars/truss3dnl.C
    Bars/truss3dnl2.C
    Bars/qtruss1d.C 
    Beams/libeam2d.C 
    Beams/libeam2dnl.C
    Beams/libeam3d.C 
    Beams/libeam3dnl.C 
    Beams/libeam3dnl2.C 
    Beams/libeam3d2.C
    Beams/beambaseelement.C
    Beams/beam2d.C 
    Beams/beam3d.C
    Beams/libeam3dboundary.C
    Beams/libeam3dboundarytruss.C
    Beams/libeam3dboundarymembrane.C
    Beams/libeam3dboundarybeam.C
    Beams/libeam3dboundaryplate.C
    Beams/libeam3dboundaryvoigt.C
    Plates/cct.C
    Plates/dkt.C 
    Plates/qdkt.C
    Plates/quad1mindlin.C
    Plates/dkt3d.C
    Shells/rershell.C 
    Shells/tr_shell01.C 
    Shells/tr_shell02.C 
    Shells/tr_shell03.C
    Shells/tr_shell11.C
    Shells/cct3d.C
    Shells/quad1mindlinshell3d.C
    Shells/shell7base.C
    Shells/tr2shell7.C
    Shells/shell7basexfem.C
    Shells/tr2shell7xfem.C
    Shells/mitc4.C
    Axisymmetry/q4axisymm.C
    Axisymmetry/axisymm3d.C 
    Axisymmetry/l4axisymm.C
    PlaneStrain/quad1planestrain.C 
    PlaneStrain/qplanestrain.C 
    PlaneStrain/qtrplanestrain.C 
    PlaneStrain/trplanestrain.C 
    PlaneStress/planestresselementevaluator.C
    PlaneStress/trplanrot.C
    PlaneStress/trplanstrss.C 
    PlaneStress/qplanstrss.C 
    PlaneStress/q9planstrss.C 
    PlaneStress/qtrplstr.C 
    PlaneStress/planstrss.C
    PlaneStress/planstrssxfem.C
    PlaneStress/trplanstrssxfem.C
    PlaneStress/qtrplanstrssxfem.C
    PlaneStress/trplanrot3d.C
    PlaneStress/trplanestressrotallman.C
    PlaneStress/basicelement.C
    PlaneStress/trplanestressrotallman3d.C
    PlaneStress/linquad3d_planestress.C
    PlaneStress/qtrplstrslip.C
    PlaneStress/qplanstrssslip.C
    AbaqusUserElement.C
    AbaqusUserElement1d.C
    AbaqusUserElement6d.C
    htselement.C
    lumpedmasselement.C
    springelement.C
    springelement3D.C
    LatticeElements/latticestructuralelement.C
    LatticeElements/lattice2d.C
    LatticeElements/lattice2dboundary.C
    LatticeElements/lattice3d.C
    LatticeElements/lattice3dboundary.C
    LatticeElements/lattice3dboundarytruss.C
    LatticeElements/latticelink3d.C
    LatticeElements/latticelink3dboundary.C
    LatticeElements/latticebeam3d.C
    LatticeElements/latticebeam3dboundary.C
    tet21ghostsolid.C
    quad1platesubsoil.C 
    quad2platesubsoil.C
    tria1platesubsoil.C
    tria2platesubsoil.C
    linedistributedspring.C
    trwarp.C
    PlaneStress/basicelement.C
    nodalspringelement.C
    #Gradient damage elements
    GradientDamage/graddamageelement.C 
    GradientDamage/Bars/truss1dgraddamage.C
    GradientDamage/Bars/qtruss1dgraddamage.C
    GradientDamage/PlaneStrain/quad1planestraingraddamage.C 
    GradientDamage/PlaneStrain/qplanestraingraddamage.C
    GradientDamage/PlaneStrain/qtrplanestraingraddamage.C
    GradientDamage/PlaneStress/trplstrgraddamage.C 
    GradientDamage/PlaneStress/qtrplstrgraddamage.C
    GradientDamage/PlaneStress/planestressgraddamage.C
    GradientDamage/PlaneStress/qplanestressgraddamage.C
    GradientDamage/3D/qspacegraddamage.C
    GradientDamage/3D/qtrspacegraddamage.C
    GradientDamage/3D/qspacegraddamage.C
    GradientDamage/3D/qwedgegraddamage.C
    #Mixed pressure elements
    MixedPressure/basemixedpressureelement.C
    MixedPressure/PlaneStrain/qtrplanestrainp1.C
    MixedPressure/PlaneStrain/quad1planestrainp1.C

    structural2delement.C
    structural3delement.C
    #phasefieldelement.C
    #qplanstrssphf.C
    )

set (sm_interface_elements
    Interfaces/cohsur3d.C
    Interfaces/structuralinterfaceelement.C
    Interfaces/intelpoint.C
    Interfaces/intelline1.C
    Interfaces/intelline1intpen.C
    Interfaces/intelline2.C
    Interfaces/intelline2intpen.C
    Interfaces/intelsurftr1.C
    Interfaces/intelsurfquad1.C
    Interfaces/bondlink3d.C
    Interfaces/bondlink3dboundary.C
    #    bondlink3dboundarybeam.C	
    #    bondlink3dboundaryplate.C		
    Interfaces/Deprecated/interfaceelement1d.C
    Interfaces/Deprecated/interfaceelem2dlin.C
    Interfaces/Deprecated/interfaceelem2dquad.C
    Interfaces/Deprecated/interfaceelem3dtrlin.C
)



set (sm_elements
    ${sm_element}
    ${sm_interface_elements}
)

add_library(sm_elements
    ${LIB_TYPE}
    ${sm_elements}
    )