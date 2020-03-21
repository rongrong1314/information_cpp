function [submap] = get_downsampled_submap(height, submap) 
    %���ݲ�ͬ�߶ȣ�������ͼ��
    [resy,resx] = size(submap);
    factor = get_resolution_from_height(height);
    %��������ڲ�ֵ�㷨�ı�ͼ��ߴ�
    tmp_map=imresize(submap,1/factor);
    submap = imresize(tmp_map,[resy,resx],'nearest');
end