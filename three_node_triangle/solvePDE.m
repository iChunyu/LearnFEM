% Solve FEM using MATLAB PDE toolbox

% XiaoCY 2019-11-30

%% main
clear;clc
close all

% create model
model = createpde('structural','static-planestrain');
gd = [2 3 ...
    0 0 6 ...
    10 0 0]';
g = decsg(gd);
geometryFromEdges(model,g);

figure
pdegplot(model,'EdgeLabels','on','VertexLabels','on')
grid on

% add properties
structuralProperties(model,'YoungsModulus',3e10,'PoissonsRatio',0.3);

% boundary condition
% structuralBC(model,'Edge',2,'XDisplacement',0);
structuralBC(model,'Edge',2,'YDisplacement',0);
structuralBC(model,'Vertex',2,'XDisplacement',0);
q = @(loc,states) 1e4*(10-loc.y);
structuralBoundaryLoad(model,'Edge',1,'Pressure',q);


% solve FEM
generateMesh(model,'GeometricOrder','linear','Hmin',2);
R = solve(model);

% results
figure
pdemesh(model)
grid on

figure
pdeplot(model,'XYData',R.Strain.exx,'ColorMap','parula')
axis equal
grid on

figure
pdeplot(model,'XYData',R.Displacement.ux,'ColorMap','parula')
axis equal
grid on