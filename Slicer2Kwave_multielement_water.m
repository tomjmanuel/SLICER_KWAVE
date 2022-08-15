% Tom 5/15/20
% script for simulating multielement arrays (not using aberration
% correction)

% xdcr files are generated by makeNiftiTransducer
% element locations can be generated by makeFibonocciArray or make
% NiftiXdcr or make IGTNiftiTransducer
clear all
%% define inputs
% output filename
fnout = 'compare_1mm';

% xdcr filenames
% point to a nifti generated by makeNiftiTransducer
% this nifti should have been used to resample CT image in slicer
xdfn = 'IGT.nii.gz';

% point to a .mat file with pixel element coordinates
elem_locs = load('IGT_pixcoords.mat');

% transducer settings
f0 = 0.65E6;               % frequency [MHz]
Amp = 1000;                % amplitude at element surface
ncyc = 5;                  % number of cycles
vox = 1E-3.*(0.25);

%% load data and set final sim size (to include padding for pml & fft efficiency)
xdcrinfo = niftiinfo(xdfn);
xdcr = niftiread(xdfn);                          
                       
% kwave needs padding around the grid for the "perfectly matched layer"
% also, kwave is much faster if you use numbers with small prime factors
dim = xdcrinfo.ImageSize;
'Your simpace size'
dim
checkFactors(min(dim),max(dim)+100)

%% look at the dimensions of your image and select a final grid size that fits
% your simspace and is in one of the first 3 rows output by checkFactors
padsize = [256 256 256];    

%% setup simulation
% put element locations into a grid
% note their vectorized order with ordervar and orderV
Egrid = zeros(dim);
ordervar = zeros(dim);
A = elem_locs.A;
for i=1:length(A)
    Egrid(A(i,1),A(i,2),A(i,3))=1;
    ordervar(A(i,1),A(i,2),A(i,3))=i;
end
orderV = ordervar(:);
orderV(orderV == 0)=[];

% pad arrays
padpre = floor((padsize-dim)/2);
padpost = ceil((padsize-dim)/2);
xdcr = padarray(xdcr,padpre,0,'pre');
xdcr = padarray(xdcr,padpost,0,'post');
Egrid = padarray(Egrid,padpre,0,'pre');
Egrid = padarray(Egrid,padpost,0,'post');
dim = padsize;

% %% test in water 
medium.density = 997.*ones(dim);
medium.sound_speed = 1480.*ones(dim);

%% create kgrid
kgrid = kWaveGrid(dim(1), vox(1), dim(2), vox(2), dim(3), vox(3));
[kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);

%% create pressure vector
source.p = createCWSignals(kgrid.t_array, f0, Amp, 0);
foo = round(1/(f0*kgrid.dt)); %samp/cyc
% zero out source after ncyc;
source.p(foo*ncyc:end)=0;

%% find focus position of xdcr from nifti file
[~,I] = max(xdcr(:));
[l, m, n]= ind2sub(dim,I);
focus_pos = [l, m, n];

% check to make sure focus_pos is inside dim
if min(dim-focus_pos)<1||min(focus_pos)<1
    error('focus position is not in the grid')
end


%% setup source
% xdcr elements are labeled as 1 in nifti
source.p_mask = xdcr==1;


%% setup sensor
sensor.mask=ones(size(kgrid.x));
sensor.record = {'p_max'};    

sensor_data = kspaceFirstOrder3DG(kgrid,medium,source,sensor);

% reshape sensor_data
pout = zeros([kgrid.Nx kgrid.Ny kgrid.Nz]);
curr=1;
for k=1:kgrid.Nz
    for j=1:kgrid.Ny
        for i=1:kgrid.Nx
            if sensor.mask(i,j,k)==1
                pout(i,j,k)=sensor_data.p_max(curr);
                curr=curr+1;
            end
        end
    end
end
% visualize
figure
subplot(131)
p=pout;
imagesc(squeeze(p(focus_pos(1),:,:)))
axis image
subplot(132)
imagesc(squeeze(p(:,focus_pos(2),:)))
axis image
subplot(133)
imagesc(p(:,:,focus_pos(3)))
axis image
% unpad and save pout to nifti with same properties as ct_simspace
% unpad
pout = pout(padpre(1)+1:end-padpost(1),padpre(2)+1:end-padpost(2),...
    padpre(3)+1:end-padpost(3));

xdcrinfo.Datatype = 'double';
%% write nifti using metadata from ctfile
%niftiwrite(pout,fnout,info,'compressed',true);
% sometimes niftiwrite modifies transforms
% you can just save the nifti with no info and pass it through
% slicer transforms to visualize in correct space
% will have to modify voxel sizes in slicer to do this
niftiwrite(pout,fnout,'compressed',true);