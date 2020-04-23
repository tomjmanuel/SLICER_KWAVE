%% getAcoustic properties
%4/16/20
%TODO: add attenuation, datacasting

function medium = getAcousticProperties(CT)
    %set water equal to air
    CT(CT<0)=0; 
    
    % take out electrodes
    CT(CT>2000)=2000;
    
    Dw = 997; Dair = 1.225; %kg/m3
    HUw = 0; HUair = -1000;

    k1 = (Dw-Dair)/(HUw-HUair);
    k0 = -Dw * HUair / (HUw-HUair);

    medium.density = k1.*CT + k0;
    %%%%%%%%%%%%%%% speed of sound
    % estimate bone indexes
    rel_CT = CT./max(CT(:)); %scale from 0 to 1 for SOS calc
    boneThresh = 1000;
    c_water = 1480;
    c_skull_max = 3100; %differs by 200 (tony uses 2900)

    bone_i = CT>boneThresh;

    medium.sound_speed = c_water.*ones(size(CT));
    medium.sound_speed(bone_i) = c_water+(c_skull_max-c_water)*rel_CT(bone_i);
end