function factor = get_resolution_from_height(altitude)
    %获得不同高度下的分辨率
    if altitude > 15
        factor = 2;
    else
        factor = 1;
    end
end