%GOPH 549 - Forward Modelling Field School
%Safian Omar Qureshi
%ID: 10086638
%
%Worked with Seismic Rockers: Sarah Reid, Sean Rhode, Tommy Boun, Andrew Ku

clear;
clc;
load data_dpdp1_10086638.mat;

%intermediate survey input variables
sp = xe(2) - xe(1); %electrode spacing       
survey_start = min(xe); %X coordinate survey start location
survey_end = max(xe); %horizontal length of survey 

%initializing loop controls
sp_factor = 1; %spacing factor is a step size to increase electrode distance as loop progresses 
j = 1;   

%this nested while loop works differently than the wenner due to how the
%potential and current electrode positions/distance have to be calculated. the inner 
%while loop calculates the different configurations of the dipole dipole
%array at horizontal survey location x0, prioritizing depth. x and z position vectors
%are calculated along with their minimum and maximum inside the loop. once all
%spacing factors are calculated for the depths, it breaks out of inner loop to outer,
%which then moves the horizontal position to the survey point and resets the spacing
%factor, recalculating over for the next survey point, prioritizing depth once
%more by changing spacing and then moving onto next horizontal point

while ((survey_start+(2*sp)+(sp_factor*sp))<=survey_end+sp_factor) 
    while (sp_factor<=10 && (survey_start+(2*sp)+(sp_factor*sp))<=survey_end+sp_factor)
        
        pot_pos = [survey_start+sp+(sp_factor*sp), survey_start+(2*sp)+(sp_factor*sp)]; %location of potential electrodes   
        curr_pos = [survey_start, survey_start+sp]; %location of current electrodes
             
        x_pos_vector(j,:) = (curr_pos(1)+(pot_pos(2)-curr_pos(1))/2);  %storing a new vector for x positions of where apparent resistivity was calculated
        z_pos_vector(j,:) = (-sp*(sp_factor+1))/2; %storing a new vector for z positions of where apparent resistivity was calculated
        x_pos_min = min(x_pos_vector);   %finding min/max of newly created x/z pseudo position vectors
        x_pos_max = max(x_pos_vector);
        z_pos_min = min(z_pos_vector);   
        z_pos_max = max(z_pos_vector);
        
        [apparent_rho_func] = resist_func(curr_pos,pot_pos); % 'rho_a' function to calculate apparent resistivity
        apparent_rho_vector(j,1) = apparent_rho_func; % storing apparent resistivity values in new vector, output variable
        
        sp_factor = sp_factor + 1; % modify spacing than loop back to the top
        j = j + 1; % updating j value to move onto next cell of vector to store data 
       
    end
    
    dx = x_pos_vector(2)-x_pos_vector(1); %finding interval in between the x, z position vectors 
    dz = dx;
    
    sp_factor = 1; %reset factor spacing back to one and do the next horizontal point over 
    survey_start = survey_start + sp; % updating x0 position to traverse the locations of electrodes horizontally
    
end

x = x_pos_min : dx : x_pos_max+sp_factor; %creating new vectors x, z to grid
z = z_pos_min : dz : z_pos_max;

[X,Z] = meshgrid(x,z); %creating grid of x, z positions of apparent resistivity locations

%output variable of interpolated generated apparent resistivity 
rho_a_interp_gen = griddata(x_pos_vector, z_pos_vector, apparent_rho_vector, X, Z); %using pseudo x, y locations of apparent resistivity, plotting it to grid x, z locations


%dipole-dipole array pseodo section plot 
figure(3);
plot(xe((1:end)),0,'rx');
hold on;
contourf(X, Z, rho_a_interp_gen);
title('Dipole-Dipole Array Pseudosection');
xlabel('x_p [m]');
ylabel('z_p [m]');
ylim([-inf,0]);


cb=colorbar;
ylabel(cb, 'p_a [ohm*m]');
colour_min = min(apparent_rho_vector);% colourbar min/max
colour_max = max(apparent_rho_vector);
caxis([colour_min, colour_max]);
prepfig;

rep_deps = [1,2,3,4,5,6,7,8,9,10];  %repeated depths on zp vector every 10 cells 

%apparent resistivity versus depth for the generated
% pseudosection and given pseudosection
figure(4);
plot(z_pos_vector(rep_deps), apparent_rho_vector(rep_deps), 'x');
hold on;
plot(zp(rep_deps), rho_a_p(rep_deps), 's');
plot(0,0);
title('Apparent Resistivity versus Depth - Dipole Dipole');
xlabel('z [m]');
ylabel('p_a [ohm*m]');
legend('Generated points', 'Given points');
prepfig;

%calculate RMS error
sum = 0;
rho_a_interp(isnan(rho_a_interp)) = 0;%changing nan values to 0 
rho_a_interp_gen(isnan(rho_a_interp_gen)) = 0;

for i = 1:length(rep_deps)
    k = (apparent_rho_vector(rep_deps(i)) - rho_a_p(rep_deps(i)))^2;
    sum = sum + k;
    
end
RMS_di = (sum/length(rep_deps))^(1/2); %output variable of script for checking error 

figure(5);
img = imread('interpretated_ref.png');
imshow(img);