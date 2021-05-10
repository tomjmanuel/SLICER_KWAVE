% make nifti transducer file 
% input makeSphericalRandom
% uses makemultibowl (kwave)
% after writing nifti, you must load it into slicer and adjust 
%   voxel sizes

clear all
load('IGTelemLocs.mat');

A = elemLocs;
% pcds (elem x (xyz))
P = [[36.4 0 62.1100]; [0 36.4 62.1100]; [-36.4 0 62.1100]; [0 -36.4 62.1100]];

% 20201207, changed elemlocs to match native coordinates (z mirrored)
%
radius_mm = 72;
diameter_mm = 6.6; %element diameter not transducer diameter
nE = 128; % n elements

usf = 4; % grid usf
padding = 30; % padding in all 3 dim
zpaddingextra = 30; % optional extra z padding

%
D = ceil(range(A))*usf;
D = D+padding;

%
fL = round(D/2);

%%

diam_pix = round(diameter_mm*usf);
if mod(diam_pix,2)==0
    diam_pix = diam_pix-1;
end

%
% shift A
A = A.*usf;
A(:,1) = A(:,1)+fL(1);
A(:,2) = A(:,2)+fL(2);
A = round(A);

% shift P
P = P.*usf;
P(:,1) = P(:,1) + fL(1);
P(:,2) = P(:,2)+fL(2);
P = round(P);

%% shift axially
% focus is axially located at 1 right now
% shift it up by the padding/2 + padding extraz
fL(3)=padding/2+zpaddingextra;
A(:,3)=A(:,3)+padding/2+zpaddingextra;
P(:,3)=P(:,3)+padding/2+zpaddingextra;

%D(3) should equal the max of the xdcr (top of cap) +padding/2
D(3)=max(A(:,3))+padding/2;


%%
%grid_size, bowl_pos, radius, diameter, focus_pos
[bowls,labels] = makeMultiBowl(D,A,radius_mm*usf,diam_pix,fL);

%% mark focus in bowls (not running this section to make nifti with targets
% make grid big enough to encomposs focus + padding
bs=3; % size of spot at focus
bowls(fL(1)-bs:fL(1)+bs,fL(2)-bs:fL(2)+bs,fL(3)-bs:fL(3)+bs)=200;
bowls(fL(1),fL(2),fL(3))=255; % mark the exact focus voxel as 255

%%
% mark a different points around focus
% this is to help create tranform that brings nifti into slicer
% using optical tracking
% assuming coordinates in bowls and labels are order xyz (so 500 steers in
% dim 1)

%called 500 but actually has points at
% 500, 090, 00-8
bowls_with500 = bowls;

%mark 50n5
dx = 1E-3/usf;
fivemmpx = round(5E-3/dx);
loc500 = fL;
loc500(1) = loc500(1)+fivemmpx;
loc500(3) = loc500(3)-fivemmpx;

bs = 2;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;
% write out bowls with 500 marked

%mark 05n5
fivemmpx = round(5E-3/dx);
loc500 = fL;
loc500(2) = loc500(2)+fivemmpx;
loc500(3) = loc500(3)-fivemmpx;
bs = 2;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;

%mark 00-5
fivemmpx = round(5E-3/dx);
loc500 = fL;
loc500(3) = loc500(3)-fivemmpx;
bs = 3;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;

% %mark 10 0n3
fivemmpx = round(5E-3/dx);
threpx = round(3E-3/dx);
loc500 = fL;
loc500(1) = loc500(1)+2*fivemmpx;
loc500(3) = loc500(3)-threpx;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;

niftiwrite(bowls_with500,'IGTwithtargs.nii','compressed',true);

%%

% outputs
% pixel coordinates to be used with timereversal
save('IGT_pixcoords.mat','A');
save('IGT_pixcoords_pcds.mat','P');
% labeled bowls
niftiwrite(labels,'IGT_labeled.nii','compressed',true);
% map with all sensor points marked 1, and focus labeled as 255, with
% surrounding voxels labeled 200
niftiwrite(bowls,'IGT.nii','compressed',true);