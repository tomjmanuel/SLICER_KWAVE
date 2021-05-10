%% This script simulates steering in water
% by computing phases and amps with Ebinni

%% load in data from Slicer2Kwave_timeReversal.m
clear all

% .mat file with element locations [mm]
elem_fn = 'manysmall_elemLocs.mat';

% .mat file with pixel coordinates [px]
pixfn = 'manysmall_pixcoords.mat';

% point at the nifti with labeled elements
fnlabeledxdcr = 'manysmall_labeled.nii.gz';

% point at nifti with unlabeled elements
% the part that emits sound should have pixels==1
fnxdcr = 'manysmall.nii.gz';

%% setup steering that you want to sim
% name for .mat that will hold pressure output
fnout= {};
fnout(1) = {'prms_ms_00n15.mat'};
fnout(2) = {'prms_ms_00n10.mat'};
fnout(3) = {'prms_ms_0010.mat'};
fnout(4) = {'prms_ms_0015.mat'};
fnout(5) = {'prms_ms_0020.mat'};
fnout(6) = {'prms_ms_00n20.mat'};
fnout(7) = {'prms_ms_500.mat'};
fnout(8) = {'prms_ms_1000.mat'};
fnout(9) = {'prms_ms_1500.mat'};
fnout(10) = {'prms_ms_2000.mat'};
fnout(11) = {'prms_ms_2500.mat'};
fnout(12) = {'prms_ms_000.mat'};

Steers = zeros([12,3]);
Steers(1,:)=[0 0 -15];
Steers(2,:)=[0 0 -10];
Steers(3,:)=[0 0 10];
Steers(4,:)=[0 0 15];
Steers(5,:)=[0 0 20];
Steers(6,:)=[0 0 -20];
Steers(7,:)=[5 0 0];
Steers(8,:)=[10 0 0];
Steers(9,:)=[15 0 0];
Steers(10,:)=[20 0 0];
Steers(11,:)=[5 0 0];
Steers(12,:)=[0 0 0];

%% loop through steering conditions
for ii=1:size(Steers,1)

    %% a few other inputs to set
    freq = 650E3;           %[Hz]
    Amp = 10000;
    steering = Steers(ii,:); %[mm] steering coords ([lat, lat, axial])
    ROC = 150;              %[mm] radius of curvature
    dx = .5;                %[mm] vox size

    %% load in element locations and compute amps and phases
    A = load(elem_fn);
    A = A.A;

    % shift axially so that focus is at zero for RS
    Ars = A;
    Ars(:,3)=Ars(:,3)-ROC;

    % compute wavenumber
    kwavenum = 2*pi*freq/1480*1e-3;

    uamp = find_uamp_RS(steering,1,kwavenum,Ars);
    Amps = abs(uamp);
    phase = angle(uamp);

    %% prepare simulation grid
    xdcr = niftiread(fnxdcr); 
    dim = size(xdcr);
    'Your simpace size'
    dim
%     checkFactors(min(dim),max(dim)+100)

    %% look at the dimensions of your image and select a final grid size that fits
    % your simspace and is in one of the first 3 rows output by checkFactors
    % make sure to give extra room if steering axially
    padsize = [540 540 432];    

    %% setup simulation
    % put element locations into a grid
    % note their vectorized order with ordervar and orderV
    Egrid = zeros(dim);
    ordervar = zeros(dim);
    A = load(pixfn); % load in pixel coordinates
    A = A.A;
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

    % test in water 
    medium.density = 997.*ones(dim);
    medium.sound_speed = 1480.*ones(dim);

    %% create kgrid
    vox = 1E-3.*[dx dx dx];
    kgrid = kWaveGrid(dim(1), vox(1), dim(2), vox(2), dim(3), vox(3));
    [kgrid.t_array, ~] = makeTime(kgrid, medium.sound_speed);

    % find focus position of xdcr from nifti file
    [~,I] = max(xdcr(:));
    [l, m, n]= ind2sub(dim,I);
    focus_pos = [l, m, n];

    % adjust focus_pos by steering amount;
    focus_pos = focus_pos + steering./dx;

    %% make a pvec for each element
    Ncyc = 5;
    t = kgrid.t_array;
    dt = kgrid.dt;
    nptspulse= round( Ncyc * 1 / (dt*freq) ); % npts in our pulse
    tp = dt.*(1:1:nptspulse); % time vector for pulse

    Ne = length(Amps);
    win = gausswin(nptspulse); %guassian window
    pvec = zeros([Ne length(t)]);

    for i=1:Ne
        pulse = Amp*Amps(i).*sin(2*pi*freq*tp+phase(i)); %pure ncycle sin wave
        pvec(i,1:nptspulse)=pulse.*win';
    end

    %% labels and mask from niftis
    labels = niftiread(fnlabeledxdcr);
    labels = padarray(labels,padpre,0,'pre');
    labels = padarray(labels,padpost,0,'post');

    % zero xdcr focus points
    xdcr(xdcr>1)=0;
    source.p_mask = xdcr;


    %% power elements with reversed signals
    % 
    labels = labels(:);
    labels(labels==0)=[];

    source.p = zeros([length(labels), size(pvec,2)]);
    for i = 1 : Ne
        % labels ==1 needs phase from sensor point where orderV ==i
        % this grabs index where orderV ==i
        [~,ind] = min((orderV-i).^2);
        % this is the phase for performing kWave based abberation correction
        source.p(labels == i, :) = repmat(pvec(ind,:), [sum(labels == i), 1]);
    end



    %% run simulation

    %sensor.mask=ones(size(kgrid.x));
    sensor.mask = zeros(size(kgrid.x));
    sensor.mask(focus_pos(1),:,:)=1;
    sensor.mask(:,focus_pos(2),:)=1;
    sensor.mask(:,:,focus_pos(3))=1;

    sensor.record = {'p_rms'};

    input_args = {'PlotSim', false, 'PMLInside', true, 'PlotPML', false, 'DataCast','single'};
    sensor_data = kspaceFirstOrder3DG(kgrid, medium, source, sensor, input_args{:});

    % reshape sensor_data
    pout = zeros([kgrid.Nx kgrid.Ny kgrid.Nz]);
    curr=1;
    for k=1:kgrid.Nz
        for j=1:kgrid.Ny
            for i=1:kgrid.Nx
                if sensor.mask(i,j,k)==1
                    pout(i,j,k)=sensor_data.p_rms(curr);
                    curr=curr+1;
                end
            end
        end
    end

    % grab planes and save them
    yz = pout(focus_pos(1),:,:);
    xy = pout(:,:,focus_pos(3));
    xz = pout(:,focus_pos(2),:);
    pdata = [];
    pdata.yz = yz;
    pdata.xy = xy;
    pdata.xz = xz;

    save(fnout{ii},'pdata');
    
end


%% functions
function uamp = find_uamp_RS(focalPoints,focalPvals,kwavenum,ucenters_mm)

pressure_vals = focalPvals;
M = length(pressure_vals);


N = length(ucenters_mm);
H = zeros([M,N]);
for m = 1:M
    Rmnvec = sqrt((ucenters_mm(:,1)-focalPoints(m,1)).^2+(ucenters_mm(:,2)-focalPoints(m,2)).^2+(ucenters_mm(:,3)-focalPoints(m,3)).^2);
    H(m,:) = 1i*exp(-1i*kwavenum*Rmnvec)./Rmnvec;
end
Hadj = conj(H.');
HHa_inv = pinv(H*Hadj);
uamp = Hadj*HHa_inv*focalPvals;
end
