% altitude = 1:0.01:26;
% var = 0.5 .* (1 - exp(-0.01 .* altitude));%������ģ��
% plot(altitude, var)
% 
%  title('Sensor model')
% xlabel('Altitude (m)')
% ylabel('Var.')
% grid minor
altitude = 1:0.01:26;
var = 0.05 .* (1 - exp(-0.2 .* altitude));
plot(altitude, var)
title('Sensor model')