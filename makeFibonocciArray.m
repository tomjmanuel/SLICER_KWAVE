% https://stackoverflow.com/questions/9600801/evenly-distributing-n-points-on-a-sphere
% https://www.sciencedirect.com/science/article/abs/pii/0025556479900804?via%3Dihub

%% create a fibonocci spiral array
R = 55;     %ROC (mm)
D = 109;    %Diam (mm)
nE = 335;   % # of elements
fnout = 'skullscan2.mat';

h = sqrt(R^2 - (D/2)^2);        % height of base of spherical cap
Acap = 2*pi*R*h;
Asphere = 4*pi*(R^2);
rat = Asphere/Acap;
npts = round(nE*rat*.95);     %npts is how many points go on the sphere, will increment up
nCap = 0;                        %holds how many points fall on spherical cap

while (nCap<nE)
    npts = npts+1;
    
    %generate points on the sphere
    [x,y,z] = getFibSphere(npts);
    x = x.*R; y = y.*R; z = z.*R;
    % crop out points that aren't on the spherical cap
    zz = z(z>h); xx = x(z>h); yy = y(z>h);
    nCap = length(zz);
    
end

if (nCap~=nE)
    'nCap ~= nE'
    nCap
end


% save
zz = -1.*zz+max(zz);
scatter3(xx,yy,zz)
axis equal
A = [xx' yy' zz'];
save(fnout,'A');

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