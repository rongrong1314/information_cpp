function [submap] = get_downsampled_submap(altitude, submap) 
    %根据不同高度，来缩放图像
    res_factor = get_resolution_from_height(altitude);
    %采用最近邻插值算法改变图像尺寸
    tmp_map = imresize(submap, 1/res_factor);
    submap = tmp_map;
end