function [submap] = get_downsampled_submap(altitude, submap) 
    %���ݲ�ͬ�߶ȣ�������ͼ��
    [res_y, res_x] = size(submap);
    factor = get_resolution_from_height(altitude);
    %��������ڲ�ֵ�㷨�ı�ͼ��ߴ�
    tmp_map = imresize(submap,1/factor);
    submap = imresize(tmp_map,[res_y, res_x],'nearest');
end