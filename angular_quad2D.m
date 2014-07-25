function [num_rays,direction_cosines,point_weights,level_weights] = angular_quad3D(N)
%Input:
% N: Order of the quadrature = # of polar angles [0,pi] 
% only defined for N <= 12
%Output:
% num_rays
% direction_cosines 
% point_weights
% level_weights
%References: B.G. Carlson 1963, Bruls 1999
%For 2D, we assume the physical domain is homogeneous in the z-axis.
%Therefore, we only do the positive hemisphere about the midplane, where
%the z-axis is the polar axis. 
%In 2D,
num_rays = N*(N+2)/2; %Up to 184
num_rays_per_octant = num_rays/4;
direction_cosines = zeros(num_rays,3);
point_weights = zeros(1,num_rays);
level_weights = zeros(1,num_rays); 

%This is the subset of [-1,1] that the direction_cosines are taken from in each
%dimension. Symmetry in each octant => only need one array
%Antisymmetric about zero
ordinates = zeros(1,N/2); 

%Choice of first cosine is arbitary. In analogy to Gaussian quadrature:
ordinates(1) = sqrt(1/(3*(N-1))); 
%In 3D space,
delta_cos =(1-3*ordinates(1)^2)* 2/(N-2); 

for i=2:N/2
    ordinates(i) = sqrt(ordinates(1)^2 + (i-1)*delta_cos);
end

%To select valid rays for the quadrature, we follow the index selection
%rule:
%sigma = i + j + k, indices of the direction cosines must equal ntheata/2 +2
%for normalization in 3D space
ray_count = 0;
for i=1:N/2
   for j=1:(N/2+1-i)
       %how do we order the rays in each octant
       ray_count = ray_count +1;
       direction_cosines(ray_count,1) = ordinates(i);
       direction_cosines(ray_count,2) = ordinates(j);
       direction_cosines(ray_count,3) = ordinates(N/2+2-i-j);
   end
end
assert(ray_count == num_rays_per_octant);

%Symmetry: reflect across second axis 
%Must be a better way of doing this
direction_cosines(ray_count+1:2*ray_count,1) = -direction_cosines(1:ray_count,1);
direction_cosines(ray_count+1:2*ray_count,2) = direction_cosines(1:ray_count,2);
direction_cosines(ray_count+1:2*ray_count,3) = direction_cosines(1:ray_count,3);
ray_count = ray_count*2;
%Symmetry: reflect across first axis
direction_cosines(ray_count+1:2*ray_count,1) = direction_cosines(1:ray_count,1);
direction_cosines(ray_count+1:2*ray_count,2) = -direction_cosines(1:ray_count,2);
direction_cosines(ray_count+1:2*ray_count,3) = direction_cosines(1:ray_count,3);
ray_count = ray_count*2;
%Material symmetry: do not reflect across final axis

%test normalization
%some of these norms have a deviation from 1.0 in the 15th decimal place
assert(all((abs(sum((direction_cosines.*direction_cosines),2) - ones(num_rays,1)) < 1e-12)));

%test output in unit sphere
quiver3(zeros(ray_count,1),zeros(ray_count,1),zeros(ray_count,1), ...
    direction_cosines(1:ray_count,1),direction_cosines(1:ray_count,2), ...
    direction_cosines(1:ray_count,3)); 

%Derive weights for integration

end