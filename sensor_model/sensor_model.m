function [var] = sensor_model(altitude, planning_parameters)
% Inverse model of sensor to detect % vegetation cover in environment.
% (Height-dependant)
%
% Input:
% altitude = current UAV altitude
% ---
% Output:
% var = variance associated with measurement

var = planning_parameters.sensor_coeff_A .* ...
    (1 - exp(-planning_parameters.sensor_coeff_B .* altitude));

 end