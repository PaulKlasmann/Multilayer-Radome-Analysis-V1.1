% Program to calculate Transmission Loss, Insertion Phase Delay and Reflection
% for a multilayer radome. 
% This program is a translation from 'Analysis of Radome-Enclosed Antennas'
% by Dennis J. Kozakoff, 1997, Artech House.
%
% Author: Dr Paul Klasmann
% Date: 05/11/2020      Version 1.1
%
% 1) Set N to equal the mumber of layers.
% 2) Set fl and fh to the lower and upper analysis frequency in GHz.
% 3) Set no_of_freq_points to a sensible value for the frequency range.
% 4) Setup the array values for the relative dielectric constant (er).
% 5) Setup the array values for the loss tangents (tand).
% 6) Setup the array values for the radome layer thicknesses (t) in mm.

close all
clear
clc
% Enter user variables here...
N = 5;                                    % Enter number of layers
fl = 17;                                  % Enter lower frequency in GHz
fh = 32;                                  % Enter upper frequency in GHz
no_of_freq_points = 301;                  % Enter No. of frequency points
no_inc_ang = 3;                           % No. of incidence angles
ang_inc_step = 15;                        % Increment angle of incidence
% Enter the dielectric contatant/permitivity for each layer
er1 = 4.0;  er2 = 1.1;  er3 = 4.0;  er4 = 1.1;  er5 = 4.0;
% Enter the loss tangent/tand for each layer
tand1 = 0.003;  tand2 = 0.001;  tand3 = 0.003;  tand4 = 0.001;  tand5 = 0.003;
% Enter the layer thickness in mm for each layer
t1 = 0.24;  t2 = 2.1;  t3 = 0.48;  t4 = 2.1;  t5 = 0.24;

er = [er1, er2, er3, er4, er5]                % Enter dielectric constants
tand = [tand1, tand2, tand3, tand4, tand5]    % Enter tand values (loss tangent)
t = [t1, t2, t3, t4, t5]                      % Enter layer thickness in mm
% End of user variables

finc = (fh-fl)/no_of_freq_points;   % Frequency increment
eo = 8.854e-12;                     % Permativity of free space (F/m)
RAD = pi/180;                       % For degree to radians conversion

% Begin Calculation
for pol = 0:1             % Pol=0 for Perpendicular Polarization
  array_index = 0;        % Pol=1 for parallel polarization
  for f = fl:finc:fh
    array_index = array_index + 1;
    LAM = 300/f;                  % Lambda is wavelength in mm
    Ko = (2*pi/LAM);              % Ko is wavenumber in mm^(-1)
    ANGLE = 0;                    % Initialize first incident angle to calculate
        for angle_index = 1:no_inc_ang;
        theta = ANGLE * RAD;      % Theta is angle of incidence in radians
        % Call function sub_wall to perform main calculations
        [Tw, IPD, Rw] = sub_wall(theta, f, er, tand, Ko, t, N, pol);
        TwdB(array_index, angle_index) = 20 * log10(Tw);
        IPDdeg(array_index, angle_index) = -IPD / RAD;
        RwdB(array_index, angle_index) = 20 * log10(Rw);
        if(IPDdeg(array_index, angle_index) == -360)
          IPDdeg(array_index, angle_index) = 0;
        end
        ANGLE = ANGLE + ang_inc_step;        % Increment angle of incidence
      end  % END ANGLE LOOP
  end      % END f LOOP

f = fl:finc:fh;  
% Transmission Loss for perpendicular polarisation
if (pol == 0)
  % Plot TL for perpendicular polarisation
  figure
  for angle_index = 1:no_inc_ang;
  plot(f, TwdB(1:end,angle_index))
  hold on;
  end
  grid on
  set(gca,'linewidth',2, 'fontsize', 14, 'XTick', fl:1:fh)
  legendtext = ['0 degrees'; '15 degrees'; '30 degrees'; '45 degrees'; '60 degrees'; '75 degrees'];
  legend( legendtext, 'location', 'southeast' );
  title('Transmission Loss for Perpendicular Polarization','FontSize', 16);
  xlabel('Frequency (GHz)','FontSize', 14);
  ylabel('Transmission Loss (dB)','FontSize', 14);
  drawnow

  % Plot IPD for perpendicular polarisation
  figure
  for angle_index = 1:no_inc_ang;
  plot(f, IPDdeg(1:end,angle_index))
  hold on;
  end
  grid on
  set(gca,'linewidth',2, 'fontsize', 14, 'XTick', fl:1:fh)
  legendtext = ['0 degrees'; '15 degrees'; '30 degrees'; '45 degrees'; '60 degrees'; '75 degrees'];
  legend( legendtext, 'location', 'southeast' );
  title('Insertion Phase Delay for Perpendicular Polarization','FontSize', 16);
  xlabel('Frequency (GHz)','FontSize', 14);
  ylabel('Insertion Phase Delay (degrees)','FontSize', 14);
  drawnow

  % Plot Reflection for perpendicular polarisation
  figure
  for angle_index = 1:no_inc_ang;
  plot(f, RwdB(1:end,angle_index))
  hold on;
  end
  grid on
  set(gca,'linewidth',2, 'fontsize', 14, 'XTick', fl:1:fh)
  legendtext = ['0 degrees'; '15 degrees'; '30 degrees'; '45 degrees'; '60 degrees'; '75 degrees'];
  legend( legendtext, 'location', 'southeast' );
  title('Reflection for Perpendicular Polarization','FontSize', 16);
  xlabel('Frequency (GHz)','FontSize', 14);
  ylabel('Return Loss (dB)','FontSize', 14);
  drawnow

% Transmission Loss for Parallel Polarisation
elseif pol == 1
  % Plot TL for parallel polarisation
  figure
  for angle_index = 1:no_inc_ang;
  plot(f, TwdB(1:end,angle_index))
  hold on;
  end
  grid on
  set(gca,'linewidth',2, 'fontsize', 14, 'XTick', fl:1:fh)
  legendtext = ['0 degrees'; '15 degrees'; '30 degrees'; '45 degrees'; '60 degrees'; '75 degrees'];
  legend( legendtext, 'location', 'southeast' );
  title('Transmission Loss for Parallel Polarization','FontSize', 16);
  xlabel('Frequency (GHz)','FontSize', 14);
  ylabel('Transmission Loss (dB)','FontSize', 14);
  drawnow

  % Plot IPD for parallel polarisation
  figure
  for angle_index = 1:no_inc_ang;
  plot(f, IPDdeg(1:end,angle_index))
  hold on;
  end
  grid on
  set(gca,'linewidth',2, 'fontsize', 14, 'XTick', fl:1:fh)
  legendtext = ['0 degrees'; '15 degrees'; '30 degrees'; '45 degrees'; '60 degrees'; '75 degrees'];
  legend( legendtext, 'location', 'southeast' );
  title('Insertion Phase Delay for Parallel Polarization','FontSize', 16);
  xlabel('Frequency (GHz)','FontSize', 14);
  ylabel('Insertion Phase Delay (degrees)','FontSize', 14);
  drawnow

  % Plot Reflection for parallel polarisation
  figure
  for angle_index = 1:no_inc_ang;
  plot(f, RwdB(1:end,angle_index))
  hold on;
  end
  grid on
  set(gca,'linewidth',2, 'fontsize', 14, 'XTick', fl:1:fh)
  legendtext = ['0 degrees'; '15 degrees'; '30 degrees'; '45 degrees'; '60 degrees'; '75 degrees'];
  legend( legendtext, 'location', 'southeast' );
  title('Reflection for Parallel Polarization','FontSize', 16);
  xlabel('Frequency (GHz)','FontSize', 14);
  ylabel('Return Loss (dB)','FontSize', 14);
  drawnow
  
end
end        % END pol LOOP