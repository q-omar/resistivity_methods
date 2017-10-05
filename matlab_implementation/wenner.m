%GOPH 549 - Forward Modelling Field School
%Safian Omar Qureshi
%ID: 10086638
%
%Worked with Seismic Rockers: Sarah Reid, Sean Rhode, Tommy Boun, Andrew Ku

clear; 
clc;
load data_wen1_10086638.mat;
      
%intermediate survey input variables
sp = xe(2) - xe(1); %electrode spacing       
survey_start = min(xe); %X coordinate survey start location
survey_end = max(xe); %horizontal length of survey 

%initializing loop controls
sp_factor = 1; %spacing factor is a step size to increase electrode distance as loop progresses
j = 1;

%this nested while loop first works traverses horizontally across the
%survey line, calculating apparent resistivities between the electrodes
%spacing at a certain depth due to a certain spacing. It also calculates the
%x and z pseudo positions and the min/max values of the
%newly created x/z position vectors
%Upon reaching the end of the survey line, the
%spacing factor is increased, the x0 counter is reset so to loop back to
%start of the line traversing horizontally once more but now at greater spacing,
%leading to greater depth. This horizontal looping resetting back to 0 with
%spacing increased keeps happening until max spacing is reached.

while (sp_factor <= 10 && (survey_start+(sp*3*sp_factor))<=survey_end+sp_factor) 
    while ((survey_start+(sp*3*sp_factor))<=survey_end+sp_factor) %nested while loop traverses horizontally through the survey
        
        pot_pos = [survey_start + (sp_factor*sp), survey_start + (sp_factor*2*sp)]; % location of potential electrodes
        curr_pos = [survey_start, survey_start + (sp*3*sp_factor)]; % location of current electrodes
        
        x_pos_vector(j,:) = (curr_pos(1)+((curr_pos(2)-curr_pos(1))/2)); %storing a new vector for x positions of where apparent resistivity was calculated
        z_pos_vector(j,:) = (-sp*sp_factor); %storing a new vector for x positions of where apparent resistivity was calculated
        x_pos_min = min(x_pos_vector); %finding min/max of newly created x/z pseudo position vectors
        x_pos_max = max(x_pos_vector);
        z_pos_min = min(z_pos_vector); 
        z_pos_max = max(z_pos_vector);
        
        [apparent_rho_func] = resist_func(curr_pos, pot_pos); % 'rho_a' function is called to calculate apparent resistivity
        apparent_rho_vector(j,1) = apparent_rho_func; % storing apparent resistivity values in new vector, output variable
        
        survey_start = survey_start + sp; % updating x0 position to traverse the locations of electrodes horizontally
        j = j + 1;    % updating j value to move onto next cell of vector to store data 
        
    end
    
    dx = x_pos_vector(2) - x_pos_vector(1); %finding interval in between the x, z position vectors 
    dz = dx;
    
    sp_factor = sp_factor + 1; %break out of nested loop once we reach end of survey to modify electrode spacing
    survey_start = 0; % recalculate survey from beginning for the updated electrode spacing 
    
end

x = x_pos_min : dx : x_pos_max+sp_factor; %creating new vectors x, z to grid
z = z_pos_min : dz : z_pos_max;

[X,Z] = meshgrid(x,z); %creating grid of x, z positions of apparent resistivity locations

%output variable of interpolated generated apparent resistivities
rho_a_interp_gen = griddata(x_pos_vector, z_pos_vector, apparent_rho_vector, X, Z); %using pseudo x, y locations of apparent resistivity, plotting it to grid x, z locations


%wenner array pseodo section plot 
figure(1);
plot(xe((1:end)),0,'rx');
hold on;
contourf(X, Z, rho_a_interp_gen);
title('Wenner Array Pseudosection');
xlabel('x_p [m]');
ylabel('z_p [m]');
ylim([-inf,0]);

cb=colorbar;
ylabel(cb, 'p_a [ohm*m]');
colour_min = min(apparent_rho_vector);% colourbar min/max
colour_max = max(apparent_rho_vector);
caxis([colour_min, colour_max]);
prepfig;

rep_deps = [1,70,136,199,259,316,370,421,469,514]; %repeated depths on zp vector

%apparent resistivity versus depth for the generated
% pseudosection and given pseudosection
figure(2);
plot(z_pos_vector(rep_deps),apparent_rho_vector(rep_deps), 'x');
hold on;
plot(zp(rep_deps), rho_a_p(rep_deps), 's');
title('Apparent Resistivity versus Depth - Wenner');
xlabel('z [m]');
ylabel('p_a [ohm*m]');
legend('Generated points', 'Given points');
prepfig;

%calculate RMS error
sum = 0;
rho_a_interp(isnan(rho_a_interp))=0; %changing nan values to 0 in given rho_a_interp data 
rho_a_interp_gen(isnan(rho_a_interp_gen))=0;

for i = 1:length(rep_deps)
    k = (apparent_rho_vector(rep_deps(i)) - rho_a_p(rep_deps(i)))^2;
    sum = sum + k;
    
end
RMS_wen = (sum/length(rep_deps))^(1/2); %output variable of script for checking error 
