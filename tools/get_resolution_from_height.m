function factor = get_resolution_from_height(height)
    %获得不同高度下的分辨率
    if height > 5
        factor = 2;
    elseif height < 1
        factor = 1;
    else
        factor = height/2;
    end
end