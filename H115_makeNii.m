% 20200503 - This script creates a nifti file for the h115 transducer,
% written by Tom. 

close all;
clear;
clc;

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
transducer_radius_m = (63.2e-03);                               % [m]
transducer_height_m = (63.2 - 54.5)*1e-03;                      % [m]
transducer_radius_pix = round(transducer_radius_m/dySim);       % [pix]
transducer_height_pix = round(transducer_height_m/dxSim);       % [pix]
sphere = makeSphericalSection(transducer_radius_pix,transducer_height_pix); 

% trasnducer has different size but same dx/dy/dz as simspace
dimTrans = size(sphere);

%% Place transducer at top of simulation space

transGrid = zeros([Nx Ny Nz]);
xoffset = 6; % [pix]
yoffset = 15;
zoffset = 15;

transGrid(xoffset:dimTrans(1)+xoffset-1,...
    yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1) = sphere;

% store focus location
xF = xoffset + transducer_radius_pix;
yF = yoffset + floor(dimTrans(2)/2);
zF = zoffset + floor(dimTrans(3)/2);

% have plane be at flat face of xdcr (will be lowest pixel after rot)
xLoc = round(xoffset + transducer_height_m/dxSim);
fLoc = [xF yF zF];

%% make grid that will be converted to nifti (makeNiftiTransducer.m)
xdcr = transGrid.*255;
bs = 3; % add points around focus to make vis easier
xdcr(fLoc(1)-bs:fLoc(1)+bs,fLoc(2)-bs:fLoc(2)+bs,fLoc(2)-bs:fLoc(2)+bs)=100;
xdcr(fLoc(1),fLoc(2),fLoc(3)) = 1; % set focus to one for future ref
xdcr(xoffset,fLoc(2),fLoc(3)) = 2; % set one pixel at top center = 2 for future ref

%% crop xdcr to match extent
xdcr = xdcr(xoffset:xoffset+xF+bs,yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1);
H115mask = xdcr;

%% write
niftiwrite(H115mask,'H115mask.nii');
