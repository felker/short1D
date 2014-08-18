function [direction_cosines,point_weights] = angular_quad1D(N)
%Input:
% N: Order of the quadrature = # of polar angles [0,pi]
%Output:
% direction_cosines 
% point_weights
%References: Davis 12
num_rays = N;
if N == 2
    direction_cosines = [1; -1];
    point_weights = [1; 1];
else
    [direction_cosines,point_weights] = lgwt(N,-1,1);
end
%test output in half circle
%quiver(zeros(N,1),zeros(N,1),direction_cosines,sqrt(ones(N,1) - direction_cosines.^2)); 

end