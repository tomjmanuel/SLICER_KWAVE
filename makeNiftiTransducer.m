% make nifti transducer file 
% input makeSphericalRandom
% uses makemultibowl (kwave)

load('transducer_32elem_R75_D75.mat')

radius_mm = 75;
diameter_mm = 6; %element diameter not transducer diameter
nE = 32; % n elements

usf = 2; % grid usf
padding = 10;

% A has element positions in mm
D = ceil(range(A))*usf;
fL = round(D/2);
D = D+padding;
fL(3)=radius_mm*usf+padding/2; % set z focus at 1;

diam_pix = diameter_mm*usf;
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
bs=2;
bowls(fL(1)+bs:fL(1)-bs,fL(2)+bs:fL(2)-bs,fL(3)+bs:fL(3)-bs)=1;

niftiwrite(labels,'trans_32elem_R75_D75_labeled.nii');

niftiwrite(bowls,'trans_32elem_R75_D75.nii');