% 20200503 - This script creates a nifti file for the h115 transducer 
% cloud of points collected using optical tracking, adapted from Tom's
% makeNii mask code. 

close all;
clear;
clc;

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

%% Setup Sim space Extent
Nx = round(max(abs(XTool))); 
Ny = round(max(abs(YTool))); 
% Nz = round(max(abs(ZTool)));

dxSim = 1e-3;
dySim = 1e-3;
dzSim = 1e-3;

% H115
transducer_radius_m = (63.2e-03);                               % [m]
transducer_height_m = (63.2 - 54.5)*1e-03;                      % [m]
transducer_radius_pix = round(transducer_radius_m/dySim);       % [pix]
transducer_height_pix = round(transducer_height_m/dxSim);       % [pix]

Nz = round(max(abs(ZTool))-transducer_radius_pix);

% Convert from mm to pixels
XToolNew = round(abs(XTool*1e-3)/dxSim) + 1; 
YToolNew = round(abs(YTool*1e-3)/dySim) + 1; 
% ZToolNew = round(abs(ZTool*1e-3)/dzSim) + 1;
ZToolNew = round(abs(ZTool+transducer_radius_pix)*1e-3/dzSim) + 1;

% Make binary mask for point cloud 
msk = zeros(max(XToolNew),max(YToolNew),max(ZToolNew));

for i = 1:length(XToolNew)
    msk(XToolNew(i),YToolNew(i),ZToolNew(i)) = 1;
end

dimTrans = size(msk);

%% Place point cloud at top of simulation space
% Not sure if should be done for point cloud. Will do it for now since we
% are also doing it for the H115 but this might change. 
transGrid = zeros([Nx Ny Nz]);
xoffset = 6; % [pix]
yoffset = 15;
zoffset = 15;

transGrid(xoffset:dimTrans(1)+xoffset-1,...
    yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1) = msk;

xdcr = transGrid.*255;

%% crop xdcr to match extent
% xdcr = xdcr(xoffset:dimTrans(1)+xoffset-1,yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1);
xdcr = xdcr(xoffset:dimTrans(1)+xoffset-1,yoffset:dimTrans(2)+yoffset-1,zoffset:dimTrans(3)+zoffset-1);
CloudPtsMask = xdcr;

%% write
niftiwrite(CloudPtsMask,'CloudPtsMask2.nii');
