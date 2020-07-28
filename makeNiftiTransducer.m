% make nifti transducer file 
% input makeSphericalRandom
% uses makemultibowl (kwave)
% after writing nifti, you must load it into slicer and adjust 
%   voxel sizes

clear all
load('george_elem_locs.mat');

radius_mm = 110;
diameter_mm = 5; %element diameter not transducer diameter
nE = 128; % n elements

usf = 4; % grid usf
padding = 30; % padding in all 3 dim
zpaddingextra = 30; % optional extra z padding

% A has element positions in mm
D = ceil(range(A))*usf;
fL = round(D/2);
D = D+padding;
fL(3)=radius_mm*usf+padding/2; % set z focus at 1;
D(3) = fL(3)+padding+zpaddingextra;
fL(1) = fL(1)+padding/2;
fL(2) = fL(2)+padding/2;

diam_pix = round(diameter_mm*usf);
if mod(diam_pix,2)==0
    diam_pix = diam_pix-1;
end

% shift A
A = A.*usf;
A(:,1) = A(:,1)+fL(1);
A(:,2) = A(:,2)+fL(2);
A(:,3) = A(:,3);
A = A+padding/2;
A = round(A);

%grid_size, bowl_pos, radius, diameter, focus_pos
[bowls,labels] = makeMultiBowl(D,A,radius_mm*usf,diam_pix,fL);

%mark focus in bowls
% make grid big enough to encomposs focus + padding
bs=5; % size of spot at focus
bowls(fL(1)-bs:fL(1)+bs,fL(2)-bs:fL(2)+bs,fL(3)-bs:fL(3)+bs)=200;
bowls(fL(1),fL(2),fL(3))=255; % mark the exact focus voxel as 255

% outputs
% pixel coordinates to be used with timereversal
save('george_pixcoords.mat','A');
% labeled bowls
niftiwrite(labels,'george_labeled.nii','compressed',true);
% map with all sensor points marked 1, and focus labeled as 255, with
% surrounding voxels labeled 200
niftiwrite(bowls,'george.nii','compressed',true);