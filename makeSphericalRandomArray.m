%https://stackoverflow.com/questions/38997302/create-random-unit-vector-inside-a-defined-conical-region/39003745#39003745
% generate points along spherical cap (random distribution)
clear all
R = 85; D = 125;
height = R-sqrt(R^2-(D/2)^2);
coneAngle  =  acos((R-height)/R);

coneDir = [0 0 1];
n=128;

z = rand(n,1)*(1-cos(coneAngle)) + cos(coneAngle);
phi = rand(n,1)*2*pi;
x = sqrt(1-z.^2).*cos(phi);
y = sqrt(1-z.^2).*sin(phi);

x = R.*x; y=R.*y; z= R.*z;

% scatter3(x,y,z)
% axis equal

% sweet now to adjust overlapping points

% go one point at a time, check its closest neighbor if its closest
% neighbor is less than a value away, try a new loc

% helps to slowly step spacing up
%
% output video
% V2 = VideoWriter('xdcrDemo');
% open(V2);
% figure 

spacing = .1:.001:5.1; %space between center of elements
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
                        tight2 = 0;
                        % get new xyz coords
                        zz = rand(1)*(1-cos(coneAngle)) + cos(coneAngle);
                        phi = rand(1)*2*pi;
                        x(i) = R*sqrt(1-zz.^2).*cos(phi);
                        y(i)= R*sqrt(1-zz.^2).*sin(phi);
                        z(i)=R*zz;
                        
                        %compare new coords to all other elements
                        for j=1:n
                            if j==i
                            else
                                dist= sqrt((x(i)-x(j)).^2+(y(i)-y(j)).^2+(z(i)-z(j)).^2);
                            
                            end
                            % if dist is smaller than spacing tolerance
                            % make a note
                            if dist<=spacing(s)
                                tight2 = 1;
                            end
                        end
                        if tight2==0
                            tight=0;
                        end
                    end
                end
            end
        end 
    end
    spacing(s)
%     if mod(s,100)==1
%         scatter3(x,y,z,'filled')
%         axis equal
%         title('Arranging elements')
%         writeVideo(V2,getframe(gcf))
%     end
end
%close(V2)

%% visualize
%close all
%[x,y,z] = sph2cart(u,v,R);
scatter3(x,y,z)
axis equal

%% save
z = -1.*(z-R);
A = [x y z];
save('transducer_128elem_R85_D125.mat','A');


