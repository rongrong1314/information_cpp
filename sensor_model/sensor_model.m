function [var] = sensor_model(altitude)
% Inverse model of sensor to detect % vegetation cover in environment.
% (Height-dependant)
%
% Input:
% altitude = current UAV altitude
% ---
% Output:
% var = variance associated with measurement

var = 0.05 .* (1 - exp(-0.2 .* altitude));

 end