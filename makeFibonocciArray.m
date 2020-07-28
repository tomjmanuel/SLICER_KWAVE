% https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
% https://www.sciencedirect.com/science/article/abs/pii/0025556479900804?via%3Dihub

%% 3d
clear all
close all
% define spherical cap
R = 110; %mm
D = 80; %mm
fn = R/D; % fnumber
h = sqrt(R^2 - (D/2)^2); % height of base of spherical cap

elemSpace = 5.56; %distance between elements [mm]

flag = 0; % flag goes to one when elements are maximally packed
npts = 3600; % n points to put on the sphere (starting number)

while ~flag

    [x,y,z] = getFibSphere(npts);
    x = x.*R; y = y.*R; z = z.*R;
    % crop out points that aren't on the spherical cap
    zz = z(z>h); xx = x(z>h); yy = y(z>h);
    scatter3(xx,yy,zz)
    axis equal
    
    space = getSpacing(xx,yy,zz);
    if space>elemSpace
        npts = npts+1;
    else
        %the maximum # of points has been found
        npts = npts-1;
        [x,y,z] = getFibSphere(npts);
        x = x.*R; y = y.*R; z = z.*R;
        % crop out points that aren't on the spherical cap
        zz = z(z>h); xx = x(z>h); yy = y(z>h);
        flag=1;
        'n elem: '
        length(zz)
        
    end
end
scatter3(xx,yy,zz)
axis equal

%% save
zz = -1.*zz+max(zz);
A = [xx' yy' zz'];
save('george_elem_locs.mat','A');

%% 2d (sunflower packing)
close all
n = 500;
indices = (0:n-1) +.5;
r = sqrt(indices/n);

theta = pi * (1 + sqrt(5)) * indices;

cosd(theta)

scatter(r.*cos(theta),r.*sin(theta))
axis equal 
%% functions to call

function [space] = getSpacing(xx,yy,zz)
    n2 = length(zz);
    dists = zeros([n2 1]);

    for i=1:n2
        dist = sqrt((xx(i)-xx).^2+(yy(i)-yy).^2+(zz(i)-zz).^2);
        dists(i) = min(dist(dist>0));
    end
    space = min(dists);
end

function [x,y,z] = getFibSphere(npts)
    n = npts;
    indices = (0:n-1) +.5;
    phi = acos(1-(2.*indices./n));
    theta = pi*(1+sqrt(5)).*indices;
    x = cos(theta).*sin(phi);
    y = sin(theta).*sin(phi);
    z = cos(phi);
end