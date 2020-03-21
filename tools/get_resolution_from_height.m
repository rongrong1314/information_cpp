function factor = get_resolution_from_height(altitude)
    %获得不同高度下的分辨率
    if altitude > 5
        factor = 2;
    elseif altitude < 1
        factor = 1;
    else
        factor = altitude/2;
    end
end