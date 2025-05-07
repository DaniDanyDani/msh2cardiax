Mesh.CharacteristicLengthMin = 0.75;
Mesh.CharacteristicLengthMax = 2.25;

Merge "epi.stl";
Merge "lv.stl";
Merge "rv.stl";
Merge "base.stl";

Coherence Mesh;

RefineMesh;

CreateTopology;


CreateGeometry;

Merge "fibrose.stl";

Surface Loop(1) = {1, 2, 3, 4, 5};

Physical Surface("base", 10) = {4};
Physical Surface("epi", 40) = {1};
Physical Surface("ve", 30) = {2};
Physical Surface("vd", 20) = {3};

Volume(1) = {1};

Surface Loop(2) = {5};

Volume(2) = {2};

Physical Volume("healthy", 1) = {1};

Physical Volume("fibrose", 2) = {2};