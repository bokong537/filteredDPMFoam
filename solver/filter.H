//- make couplingFilter class
Time coupleDiffusionRunTime
(
    "controlDiffDict",
    args
);

fvMesh coupleDiffusionMesh
(
    IOobject
    (
        fvMesh::defaultRegion,
        coupleDiffusionRunTime.timeName(),
        coupleDiffusionRunTime,
        IOobject::MUST_READ
    )
);

simpleControl coupleDiffusionSimple(coupleDiffusionMesh);

couplingFilter filterModel
(
    mesh,
    coupleDiffusionRunTime,
    coupleDiffusionMesh,
    coupleDiffusionSimple
);



