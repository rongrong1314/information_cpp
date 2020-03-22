function grid_map = ...
    update_map_with_correlation(pos, submap, ...
    grid_map, submap_coordinates, planning_parameters)
% Updates grid map at a UAV position using measurements
% received using height-dependent sensor model.
%% Update occupancy grid map - weeds.
% Compute variance of measurements.
%测量误差作为方差，老版本：草：传感器模型，非草：0.001
var = ones(size(submap))*sensor_model(pos(3), planning_parameters);
R = diag(reshape(var,[],1));
% Create measurement model.
H = construct_H(grid_map.m, submap, submap_coordinates, pos(3));

%% Correlated Fusion
% Obtain maximum aposteriori estimate using Bayesian Fusion.
x = reshape(grid_map.m,[],1);
z = reshape(submap,[],1);
P = grid_map.P;

%x为观测矩阵，z为观测值，H为状态观测矩阵，是一个单位矩阵，z-H*x为测量余量
%根据观测值预测
[x,Pf] = KF_update_cholesky(x,P,z-H*x,R,H);
grid_map.m = reshape(x, size(grid_map.m,1), size(grid_map.m,2));
grid_map.P = Pf;

end
%%
% [w_win_j, w_win_i] = find(submap == 1);
% 
% w_ind = sub2ind(size(submap), w_win_j, w_win_i);
% %子图栅格坐标
% [sub_submapx, sub_submapy] = meshgrid([submap_coordinates.xl:submap_coordinates.xr],...
%         [submap_coordinates.yd:submap_coordinates.yu]);
% %子图全局坐标
% w_ind_global = sub2ind(size(grid_map.m),reshape(sub_submapy,[],1),reshape(sub_submapx,[],1));
%   
% var = ones(size(submap));
% 
% if (~isempty(w_win_i))
%     weeds = 1;
%     %转换杂草坐标为全局坐标
%     w_x = submap_coordinates.xl + w_win_i - 1;
%     w_y = submap_coordinates.yd + w_win_j - 1;
%     w_win_ind = sub2ind(size(grid_map.m), w_y, w_x);%转换为1维
%     
%     %更新杂草位置的概率，传感器误差
%     submap(w_ind) = classify_weeds(pos_env(3), size(w_win_ind,1), planning_parameters,weeds);
%     varw = (pos_env(3)/64);
%     var(w_ind)=varw;
% end
% 
% %% Update occupancy grid map - nonweeds.
% % Find nonweeds, extract indices, and update probabilities.
% [nw_win_j, nw_win_i] = find(submap == 0);
% 
% nw_ind = sub2ind(size(submap), nw_win_j, nw_win_i);
% 
% if (~isempty(nw_win_i))
%     weeds = 0;
%     nw_x = submap_coordinates.xl + nw_win_i - 1;
%     nw_y = submap_coordinates.yd + nw_win_j - 1;
%     nw_win_ind = sub2ind(size(grid_map.m), nw_y, nw_x);
%     
%     submap(nw_ind) = classify_weeds(pos_env(3), size(nw_win_ind,1), planning_parameters,weeds);
%     varnw = 0.001;
%     var(nw_ind)=varnw;
% end
% 
% %% Perform correlated fusion.
% 
% [grid_map.m, grid_map.P] = fuse_measurements(grid_map.m, grid_map.P, ...
%     submap,var, w_ind_global);


