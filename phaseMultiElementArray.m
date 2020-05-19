%% load in data from Slicer2Kwave_timeReversal.m
clear all
load('patsesnsorwater.mat')

% point at the nifti with labeled elements
fnlabeledxdcr = 'xdcr_128elem_R100_D100_vox250um_labeled.nii.gz';

% point at nifti with unlabeled elements
% the part that emits sound should have pixels==1
fnxdcr = 'xdcr_128elem_R100_D100_vox250um.nii.gz';

% name for nifti that will hold pressure output
fnout = 'prms_testmultielement_water';

%% visualize pdata 
% use this plot to select a time window that incorporates
% the bulk of the signal arriving at the sensor
% we don't want the phase of reflected sounds later in the sim
% set t0 and te below
plot(sensor_data.p(10:10:end,:)')


%% a few inputs to set
freq = f0; %0.5E6; %[MHz] should load in with other sensor_data
t0 = 1200; % beginning time to calc phase
te = 1500; % end time to calc phase
Amp = 100;

%%
% Estimate phase for every pixel on the sensor mask
t =kgrid.t_array;
dt = kgrid.dt;
pp = sensor_data.p;
p = flip(pp(:,t0:te),2); % this is equivalent to conjugating phase
dimP=  size(p);

%demod to get I and Q (thanks for teaching me Byram)
Re = sum(pp.* repmat(cos(2*pi*freq.*t),[dimP(1) 1]),2);
Im = sum(pp.* repmat(sin(2*pi*freq.*t),[dimP(1) 1]),2);

foo = complex(Re,Im);
phase = angle(foo);
Amps = abs(foo);

%% make a pvec for each element
Ncyc = 10;
t = kgrid.t_array;
dt = kgrid.dt;
nptspulse= round( Ncyc * 1 / (dt*freq) ); % npts in our pulse
tp = dt.*(1:1:nptspulse); % time vector for pulse

Ne = length(Amps);
win = gausswin(nptspulse); %guassian window
pvec = zeros([Ne length(t)]);



for i=1:Ne


    pulse = Amp.*sin(2*pi*freq*tp+ phase(i)); %pure ncycle sin wave
    pvec(i,1:nptspulse)=pulse.*win';

   
end

%% get mask and labels from niftis
%[source.p_mask, labels] = makeMultiBowl(dim, [xx, yy, zz], radius, diameter, [focus_pos(1) focus_pos(2) focus_pos(3)], 'RemoveOverlap',true);
labels = niftiread(fnlabeledxdcr);
xdcr = niftiread(fnxdcr);

% zero xdcr focus points
xdcr(xdcr>1)=0;

% pad both grids to match sim space size
xdcr = padarray(xdcr,padpre,0,'pre');
xdcr = padarray(xdcr,padpost,0,'post');
labels = padarray(labels,padpre,0,'pre');
labels = padarray(labels,padpost,0,'post');

source.p_mask = xdcr;

%% visualize setup
% figure;
% pfoo = imresize3(source.p_mask,.2);
% cfoo = imresize3(medium.sound_speed,.2);
% 
% pfoo = flip(pfoo,3);
% cfoo = flip(cfoo,3);
% 
% isosurface(pfoo);axis equal;hold on;
% isosurface(cfoo,1600);hold on;
% 

%% power elements with reversed signals
% 
labels = labels(:);
labels(labels==0)=[];

source.p = zeros([length(labels), size(pvec,2)]);

for i = 1 : Ne
    % labels ==1 needs phase from sensor point where orderV ==1
    % this grabs index where orderV ==i
    [~,ind] = min((orderV-i).^2);
    
    source.p(labels == i, :) = repmat(pvec(ind,:), [sum(labels == i), 1]);
end



%% run simulation

sensor.mask=ones(size(kgrid.x));
% sensor.mask = zeros(size(kgrid.x));
% sensor.mask(focus_pos(1),:,:)=1;
% sensor.mask(:,focus_pos(2),:)=1;
% sensor.mask(:,:,focus_pos(3))=1;

sensor.record = {'p_max'};

input_args = {'PlotSim', false, 'PMLInside', true, 'PlotPML', false};

sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

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

info = niftiinfo(fnxdcr);
% if file won't write, you may need to change some properties in info
%info.Datatype = 'single';
%% write nifti using metadata from ctfile
niftiwrite(pout,fnout,info,'compressed',true);