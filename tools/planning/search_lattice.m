function path = search_lattice(point_init, lattice, grid_map, map_parameters, ...
    planning_parameters)
%执行贪婪网格搜索候选人名单，基于一个信息性的目标以确定最有希望的点访问。
%起点固定(此处不测量)
% ---
% Inputs:
%起点，要评估的候选，当前地图
% point_init: starting location
% lattice: list of candidates to evaluate
% grid_map: current grid map (mean + covariance)
% ---
% Output:
% path: grid search result

P_trace_prev = trace(grid_map.P);
point_prev = point_init;
path = point_init;

while (planning_parameters.control_points > size(path, 1))
    
    % Initialise best solution so far.
    obj_min = Inf;
    point_best = -Inf;
    
    for i = 1:size(lattice, 1)
        
        point_eval = lattice(i, :);
        grid_map_eval = predict_map_update(point_eval, grid_map, ...
            map_parameters, planning_parameters);
        P_trace = trace(grid_map_eval.P);
        
        gain = P_trace_prev - P_trace;
        cost = pdist([point_prev; point_eval]);
        obj = -gain*exp(-planning_parameters.lambda*cost);
        
        % Update best solution.
        if (obj < obj_min)
            obj_min = obj;
            point_best = point_eval;
        end
        
    end
    
    % Update the map with measurement at best point.
    grid_map = predict_map_update(point_best, grid_map, ...
        map_parameters, planning_parameters);
    disp(['Point ', num2str(size(path,1)+1), ' at: ', num2str(point_best)]);
    disp(['Trace of P: ', num2str(trace(grid_map.P))]);
    path = [path; point_best];
    
    P_trace_prev = trace(grid_map.P);
    point_prev = point_best;
    
end

end