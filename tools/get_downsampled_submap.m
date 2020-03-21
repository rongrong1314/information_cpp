function [submap] = get_downsampled_submap(height, submap) 
    %根据不同高度，来缩放图像
    [resy,resx] = size(submap);
    factor = get_resolution_from_height(height);
    %采用最近邻插值算法改变图像尺寸
    tmp_map=imresize(submap,1/factor);
    submap = imresize(tmp_map,[resy,resx],'nearest');
end