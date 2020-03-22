function t = estimate_segment_time(pos1, pos2, v_max, a_max, constant) 
if (nargin < 5)
  % Magic nfabian constant. :(
  constant = 6.5;
end

distance = norm(pos2 - pos1);
t = 2*distance/v_max*(1+constant*v_max/a_max * exp(-distance/v_max*2));

if (t < 0.01)
  t = 0.01;
end
end