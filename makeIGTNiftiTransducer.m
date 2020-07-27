% make nifti transducer file 
% input makeSphericalRandom
% uses makemultibowl (kwave)
% after writing nifti, you must load it into slicer and adjust 
%   voxel sizes

clear all
load('IGTelemLocs.mat');

A = elemLocs;
%
radius_mm = 72;
diameter_mm = 6.6; %element diameter not transducer diameter
nE = 128; % n elements

usf = 2; % grid usf
padding = 30; % padding in all 3 dim
zpaddingextra = 30; % optional extra z padding

%
% A has element positions in mm
D = ceil(range(A))*usf;
D = D+padding;

%
fL = round(D/2);
fL(3)=radius_mm*usf+padding/2; % set z focus at 1;
D(3) = fL(3)+padding+zpaddingextra;

diam_pix = round(diameter_mm*usf);
if mod(diam_pix,2)==0
    diam_pix = diam_pix-1;
end

%
% shift A
A = A.*usf;
A(:,1) = A(:,1)+fL(1);
A(:,2) = A(:,2)+fL(2);
A(:,3) = A(:,3)+padding/2;
%A = A+padding/2;
A = round(A);


%%
%grid_size, bowl_pos, radius, diameter, focus_pos
[bowls,labels] = makeMultiBowl(D,A,radius_mm*usf,diam_pix,fL);

%% mark focus in bowls (not running this section to make nifti with targets
% make grid big enough to encomposs focus + padding
bs=5; % size of spot at focus
bowls(fL(1)-bs:fL(1)+bs,fL(2)-bs:fL(2)+bs,fL(3)-bs:fL(3)+bs)=200;
bowls(fL(1),fL(2),fL(3))=255; % mark the exact focus voxel as 255

%%
% mark a different points around focus
% this is to help create tranform that brings nifti into slicer
% using optical tracking
% assuming coordinates in bowls and labels are order xyz (so 500 steers in
% dim 1)

% zero everythign to make rendering easier
bowls_with500 = bowls.*0;

%mark 500
dx = 1E-3/usf;
fivemmpx = round(5E-3/dx);
loc500 = fL;
loc500(1) = loc500(1)+fivemmpx;

bs = 2;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;
% write out bowls with 500 marked

%mark 090
fivemmpx = round(9E-3/dx);
loc500 = fL;
loc500(2) = loc500(2)+fivemmpx;
bs = 2;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;

%mark 00-8
fivemmpx = round(8E-3/dx);
loc500 = fL;
loc500(3) = loc500(3)-fivemmpx;
bs = 3;
bowls_with500(loc500(1)-bs:loc500(1)+bs,loc500(2)-bs:loc500(2)+bs,loc500(3)-bs:loc500(3)+bs)=200;

%mark 000
bs=1; % size of spot at focus
bowls_with500(fL(1)-bs:fL(1)+bs,fL(2)-bs:fL(2)+bs,fL(3)-bs:fL(3)+bs)=200;
bowls_with500(fL(1),fL(2),fL(3))=255; % mark the exact focus voxel as 255


niftiwrite(bowls_with500,'IGTwithtargs.nii','compressed',true);

%%

% outputs
% pixel coordinates to be used with timereversal
save('IGT_pixcoords.mat','A');
% labeled bowls
niftiwrite(labels,'IGT_labeled.nii','compressed',true);
% map with all sensor points marked 1, and focus labeled as 255, with
% surrounding voxels labeled 200
niftiwrite(bowls,'IGT.nii','compressed',true);