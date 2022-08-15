%% getAcoustic properties
%7/29/20
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
    boneThresh = 600;
    c_water = 1480;
    c_skull_max = 3100; %differs by 200 (tony uses 2900)

    bone_i = CT>boneThresh;

    medium.sound_speed = c_water.*ones(size(CT));
    medium.sound_speed(bone_i) = c_water+(c_skull_max-c_water)*rel_CT(bone_i);

        %%%%%%%%%%%%%%%% absorption
    % from manual: alpha_coeff units: [dB/Mhz/cm] (page 32)
    %              alpha_power: scales alpha with frequency (just set to 1) we
    %                   are close to 1Mhz anyway
    %                   can't actually set to 1 but set close to 1
    % Ref 2 explains this well
    % ref 1, Pinton paper about absorption, scattering, etc
    %from ref 1: bone gives absoprtion of 2.7 dB/Mhz/cm longitudinal 
    %            and 5.4 dB/Mhz/cm shear...
    % i will create binary mask of bone and set alpha to 8 dB/Mhz/cm within
    % mask, 0.2 dB/Mhz/cm everywhere else
    
    % ref 2 duck book table 4.15
    % which cites this paper 
    % http://www.brl.uiuc.edu/Publications/1978/Fry-JASA-1576-1978.pdf

    medium.alpha_power = 1.1; % from aubry and commonly used
    medium.alpha_coeff = 0.2.*ones(size(CT));
    medium.alpha_coeff(bone_i) = 2.7+5.4; % ref 1
    medium.alpha_coeff(bone_i) = 22; % ref 2 
end