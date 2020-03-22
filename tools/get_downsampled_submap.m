function [submap] = get_downsampled_submap(altitude, submap) 
    %���ݲ�ͬ�߶ȣ�������ͼ��
    res_factor = get_resolution_from_height(altitude);
    %��������ڲ�ֵ�㷨�ı�ͼ��ߴ�
    tmp_map = imresize(submap, 1/res_factor);
    submap = tmp_map;
end