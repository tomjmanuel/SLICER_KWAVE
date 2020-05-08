% 20200504 - Register cloud of points to model of H115 using the ICP
% algorithm.

close all;
clear;
clc;

addpath('./ICP_scripts');

%% Setup Sim space Extent
Nx = 80; Ny = 80; Nz = 80;

% using h115, focusing down x dim
XrangeSim = .08;
YrangeSim = .08;
ZrangeSim = .08;

dxSim = XrangeSim/Nx;
dySim = YrangeSim/Ny;
dzSim = ZrangeSim/Nz;

%% Setup transducer
% H115
transducer_radius_m = (63.2e-03);
transducer_height_m = (63.2 - 54.5)*1e-03;
transducer_radius_pix = round(transducer_radius_m/dySim);
transducer_height_pix = round(transducer_height_m/dxSim);
sphere = makeSphericalSection(transducer_radius_pix,transducer_height_pix);

% trasnducer has different size but same dx/dy/dz as simspace
dimTrans = size(sphere);

%% Read in the tsv file
% Combining two datasets so we have more points for the cloud.
D = dir('*.tsv');
XTool = [];
YTool = [];
ZTool = [];

for k = 1:length(D)
    try
        fileName = D(k).name;
        fid = fopen(fileName);
        C = textscan(fid, '%f %s %f %f %s %f %f %f %f %f %f %f %f %f %s %f %f %s %f %f %f %f %f %f %f %f %f', 'Delimiter','\t');

        foo = strsplit(C{1,2}{1,1},' ');
        tool1Name = foo{4};
        foobar = strsplit(C{1,15}{1,1},' ');
        tool2Name = foobar{4};

        XTool = [XTool; (C{1,10})];
        YTool = [YTool; (C{1,11})];
        ZTool = [ZTool; (C{1,12})];
    end
end

XToolNew = XTool - min(XTool); 
YToolNew = YTool - min(YTool); 
ZToolNew = ZTool - min(ZTool); 

plot3(XTool,YTool,ZTool,'.');  
axis('image');

%% Register using ICP Algorithm
% Find X,Y,Z coordinates for transducer
ind = find(sphere~=0);
[targetX,targetY,targetZ] = ind2sub(size(sphere), ind);

target = [targetX,targetY,targetZ];
source = [XTool YTool ZTool];

[error,Reallignedsource,transform] = rigidICP(source,target);

figure;
plot3(Reallignedsource(:,1),Reallignedsource(:,2),Reallignedsource(:,3),'.','MarkerSize',8)
hold on;
plot3(source(:,1),source(:,2),source(:,3),'.','MarkerSize',12);
axis('image');
legend('MakeSphericalSection','SurfaceData');

%% Write to text file
addpath('./Transform_scripts');

param = [transform.T(1,:) transform.c(1,1); transform.T(2,:) ...
     transform.c(1,2); transform.T(3,:)  transform.c(1,3); 0 0 0 1];
param = param.*[1 1 1 -1; -1 -1 -1 -1; -1 -1 -1 -1; 0 0 0 1];
 
% I still need to further understand how Slicer and MATLAB define their
% coordinate systems, there needs to be a 180 deg rotation to get
% everything in the right space. 
 
cli_lineartransformwrite('H115_to_CloudPt.txt', param, 'slicer')

%% 
% ninfo = niftiinfo('H115maskTest.nii');
% ninfo.Transform.T = param;
