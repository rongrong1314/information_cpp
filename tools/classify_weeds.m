function [prob] = classify_weeds(altitude, num, planning_parameters,weeds)
% Classifier sensor model for detecting weeds:
% Outputs probability of finding a weed at a given point, for current
% altitude.

if (altitude > 30)
    prob = 0.5;
else
    if(weeds == 1)
        %���߶ȴ���������ģ�ͣ����Ӳݸ���
        prob = polyval(planning_parameters.weed_coeffs, altitude);
        prob = prob*ones(num,1);
    else
       prob = polyval(planning_parameters.nonweed_coeffs, altitude);
       prob = prob*ones(num,1);
end

end

