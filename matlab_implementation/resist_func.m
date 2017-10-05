%GOPH 549 - Forward Modelling Field School
%Safian Omar Qureshi
%ID: 10086638
%
%Worked with Seismic Rockers: Sarah Reid, Sean Rhode, Tommy Boun, Andrew Ku


function [apparent_rho_func] = resist_func(pot_pos,curr_pos) %this function is similar to the one developed on field for ERT, passing potential/current electrode positions from main script

    %initial input variables are earth model parameters that are iteratively
    %tweaked to create stronger model. these input variables are resistivities rho1,rho2 of
    %surface and basal halfspace layers, interface depth z, current and
    %potential electrode locations based on given input vector xe
    
    %current is cancelled in the mathematical equation due to assumptions and not needed
    %specifically in this implementation as a variable 
    
    rho1 = 332; %first layer resistivity (ohm m) 
    rho2 = 865; %second layer resistivity
    z = 23; %interface depth (m)

    %intermediate inputs 
    AM = abs((pot_pos(1))-(curr_pos(1))); %current/potential positions parameters which are passed to function from runme script to calculate distance between them
    BM = abs((pot_pos(2))-(curr_pos(1)));
    AN = abs((pot_pos(1))-(curr_pos(2))); %and the distance between them are calculated inside this function to use for calculating apparant resisitivity values 
    BN = abs((pot_pos(2))-(curr_pos(2)));

    ref_coef = (rho2-rho1)/(rho2+rho1); % reflection coefficient
    k = 2*pi*(1/(((1/AM)-(1/AN))-((1/BM)-(1/BN)))); % geometric factor
    s = (rho1)/(k); %rearrange formula by Telford, p551 (Eq 8.37) for numerator 

    m = 1; % starting with term 1, m = 1  
    limit = 1; % limit initialized at 1
    denom_coef = 4*(m^2)*(z^2); 

    %this loop is similar to the one developed for ERT where numerator and denominator are
    %continuously calculated until a threshold is reached, at which point the loop is broken
    %to give us parameter s which calculates apparent resisitivity when multiplied by k
    
    threshold = 1e-3; %small number to compare change in ratio
    while (limit > threshold)

        denom_coef = 4*(m^2)*(z^2);
        dev_v1 = ((rho1)/(2*pi))*(2*ref_coef^m*((1/(AM^2+denom_coef)^(1/2)-1/(AN^2+denom_coef)^(1/2))-(1/(BM^2+denom_coef)^(1/2)-1/(BN^2+denom_coef)^(1/2))));  %Telford, p551 (Eq 8.37)
        s = s + dev_v1;
        limit = abs(dev_v1/s);
        m=m+1;

    end

    apparent_rho_func = s*k; %output variable of the function: apparent resisitivty calculated at a specific point 

end

    