function factor = get_resolution_from_height(altitude)
    %��ò�ͬ�߶��µķֱ���
    if altitude > 15
        factor = 2;
    else
        factor = 1;
    end
end