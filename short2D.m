%2D short characteristics solver for RT
%Time independent problem
%Main ref: Davis, et al 2012
%Other refs: Kunasz, Auer 1988, Bruls 1999, Mihalas,Auer,Mihalas 78
%Vectors are column vectors

%Parameters
clear all;
close all;
nx = 100;
ny = 100;
ntheta = 12; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lx = 1.0;
ly = 1.0;
c = 1.0;
dx = lx/nx;
dy = ly/ny;

%Angular Discretization, Sn discrete ordinates. 
[na, mu, pw, lw] = angular_quad2D(ntheta);

%Spatial Discretization
%Radiation points are centered on the fluid cells
xx=linspace(0,lx,nx)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nx,ny,na); 
mean_intensity = zeros(nx,ny); 

%Total opacity 
X_tot = zeros(nx,ny); 
%Photon destruction probability eps = X_abs/X_tot
%This is set to 1, as absorption and emission dominate in LTE
destruction_probability = ones(nx,ny); 
%Thermal blackbody source function
T = 1000; % in K
temperatures = T*ones(nx,ny); %uniform domain temperature
%Isotropic blackbody function
%B_v(T) = 2*h*v^3/c^2 * 1/(e^(hv/kT) -1)
B = 1e-3;
thermal_source = B*temperatures;

%Calculate source function
%Integrate I for mean intensity, J using quadrature rule
%mean_intensity = 1/(4*pi)*sum(weights*ones(1,ntheta)*intensity);
source_function = destruction_probability.*thermal_source; % + ...
%(ones(nx,ny) - destruction_probability).*mean_intensity;

%Optical depth intervals before and after pt: depends on opacticy and location
opt_depth = zeros(2,1);
%Interpolation coefficients
interp_coeff = zeros(3,1);

%Test Problems:
%Coherent searchlight beam from LHS of domain
%intensity(1,:) = 1;
%X_tot = zeros(nx,ny); %this implies no blackbody emission or absorption

%Uniform material
%X_tot = 10*ones(nx,ny);

for j=1:ntheta
    %Need to sweep domain angles in an organized manner.
    %We fist start in bottom LHS corner of domain, compute upward rays
    %(mu(2) > 0) and right pointing rays (mu(1) > 0) for this row
    %Then, compute left pointing rays for this row (mu(1) < 0)
    
    %Shallow rays must be handled differently for qudratic and higher
    %interpolation
    
    if mu(j) >=0
        first = 2;
        last = nz-1;
        upwind = -1;
        downwind = 1;
    else
        first = nz-1;
        last = 2;
        upwind = 1;
        downwind = -1;
    end
    for k=first:downwind:last %trace characteristics downwind 
        %Hayek10 computes the optical depth intervals using Bezier interp
        %Davis12 simply takes a linear interpolation
        %path integral from previous point 
        opt_depth(1) = dz*(X_tot(k)+X_tot(k+upwind)/2);
        %path integral to next point
        opt_depth(2) = dz*(X_tot(k+downwind)+X_tot(k)/2);
        %Calculate interpolation coefficients
        %Possible schemes: linear, parabolic, Bezier with control points
        interp_coeff = sc_interpolation(opt_depth(1),opt_depth(2));
        %Update
        intensity(k,j) = intensity(k+upwind,j)*exp(-opt_depth(1)) + interp_coeff(1)*source_function(k+upwind) + ...
            +interp_coeff(2)*source_function(k) + interp_coeff(3)*source_function(k+downwind);
     end
end

%Output
%Plot each angular intensity (recall, ntheta must be even)
for i=1:ntheta
    hi = subplot(4,3,i); 
    plot(hi,zz(2:nz-1),intensity(2:nz-1,i));
    x_label = sprintf('mu = %f',mu(i));
    xlabel(x_label);
end
    h = figure;
    quiver(zeros(ntheta,1),zeros(ntheta,1),mu,sqrt(ones(ntheta,1) - mu.^2)); 
