%https://stackoverflow.com/questions/38997302/create-random-unit-vector-inside-a-defined-conical-region/39003745#39003745
% generate points along spherical cap (random distribution)
clear all
R = 75; D = 75;
height = R-sqrt(R^2-(D/2)^2);
coneAngle  =  acos((R-height)/R);

coneDir = [0 0 1];
n=32;

z = rand(n,1)*(1-cos(coneAngle)) + cos(coneAngle);
phi = rand(n,1)*2*pi;
x = sqrt(1-z.^2).*cos(phi);
y = sqrt(1-z.^2).*sin(phi);

x = R.*x; y=R.*y; z= R.*z;

scatter3(x,y,z)
axis equal

%% sweet now to adjust overlapping points

% go one point at a time, check its closest neighbor if its closest
% neighbor is less than a value away, try a new loc

% helps to slowly step spacing up
%
spacing = .1:.0005:6.1; %space between center of elements
for s=1:length(spacing)
    for i=1:n
        for j=1:n
            if i==j
                %don't compare to itself
            else
                dist = sqrt((x(i)-x(j)).^2+(y(i)-y(j)).^2+(z(i)-z(j)).^2);
                if dist<=spacing(s)
                    tight = 1;
                    while(tight)
                        % get new xyz coords
                        zz = rand(1)*(1-cos(coneAngle)) + cos(coneAngle);
                        phi = rand(1)*2*pi;
                        x(i) = R*sqrt(1-zz.^2).*cos(phi);
                        y(i)= R*sqrt(1-zz.^2).*sin(phi);
                        z(i)=R*zz;
                        
%                         uu = rand(); vv = .125*NA*rand();
%                         thet = 2*pi*uu; phi = acos(2*vv-1);
%                         x(i) = R.*sin(phi).*cos(thet);
%                         y(i) = R.*sin(phi).*sin(thet);
%                         z(i) = R.*cos(phi);
                        %check again
                        for j=1:n
                            if j==i
                            else
                                dist= sqrt((x(i)-x(j)).^2+(y(i)-y(j)).^2+(z(i)-z(j)).^2);
                            end
                            if dist<=spacing(s)
                                break
                            end
                        end
                        if j==n
                            tight=0;
                        end
                    end
                end
            end
        end 
    end
    spacing(s)
end

%% visualize
%close all
%[x,y,z] = sph2cart(u,v,R);
scatter3(x,y,z)
axis equal

%% save
z = -1.*(z-R);
A = [x y z];
save('transducer_32elem_R75_D75.mat','A');


