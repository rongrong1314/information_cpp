function factor = get_resolution_from_height(altitude)
    %��ò�ͬ�߶��µķֱ���
    if altitude > 5
        factor = 2;
    elseif altitude < 1
        factor = 1;
    else
        factor = altitude/2;
    end
end