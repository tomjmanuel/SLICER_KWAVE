%% Use HAS to compute ideal phases for xdcr

% Inputs %%%%
% xdcr_filename
% 
% ct data

% Outputs %%%%
% [Nelem 1] phase array


% Inputs %%%%%%%%%%%%%%
xdfn = 'IGT.nii.gz';            % xdcr nifti file (pix 255 at focus)
ctfn = 'ct_simspace.nii.gz';    % simspace CT 
elemfn = 'IGT_pixcoords.mat';   % element pixel locations
f = 0.65E6;                        % frequency
%%%%%%%%%%%%%%%%%%%%%%

xdcr = niftiread(xdfn);
ct = niftiread(ctfn);
info = niftiinfo(ctfn);
load(elemfn) % should load as A

%% get acoustic maps
ct = double(ct);
if min(ct(:)>=0)
    ct=ct-1024;
    'shifting ct data'
end
medium = getAcousticProperties(ct);

%% identify focus position in xdcr
[~,ind] = max(xdcr(:));
[l,m,n] = ind2sub(size(xdcr),ind);
focus_pos = [l m n];

%% call HAS
%create virtual point source
dim = size(xdcr); 
p0 = zeros([dim(1) dim(2)]);
p0(l,m)= 100;
vox = info.PixelDimensions(1) * 1E-3;
p = zeros(dim);
%[p, Refp, Forp] = HAS(p0,vox,f,sos,density)
%[~,~,p] = HAS(p0, vox,f, flip(medium.sound_speed(:,:,1:n),3), flip(medium.density(:,:,1:n),3));
[~,~,p(:,:,n:end)] = HAS(p0, vox,f, medium.sound_speed(:,:,n:end), medium.density(:,:,n:end));
% flip dim so that HAS propogates in correct direction
%p = flip(p,3);
%% get phases
% create a 3d grid that only has phase points at the element locations
nE = size(A,1);

phasegrid = 10.*ones(dim); % set == 10 so that we can delete all 10's later

for i=1:nE
    foo = angle(p(A(i,1),A(i,2),A(i,3)));
    phasegrid(A(i,1),A(i,2),A(i,3)) = foo;
end

% vectorize and remove 0s
phaseHAS = phasegrid(:);
phaseHAS( phaseHAS==10)= [];
