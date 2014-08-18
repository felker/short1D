%1D short characteristics solver for RT
%Time independent problem
%Main ref: Davis, et al 2012
%Other refs: Kunasz, Auer 1988, Bruls 1999, Mihalas,Auer,Mihalas 78
%Vectors are column vectors

%Parameters
clear all;
close all;
nz = 100;
ntheta = 2; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].

lz = 1.0;
c = 1.0;
dz = lz/nz;

%Angular Discretization, Sn discrete ordinates. 
[mu, w] = angular_quad1D(ntheta);

%Spatial Discretization
%Radiation points are centered on the fluid cells
zz=linspace(0,lz,nz)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nz,ntheta); 
mean_intensity = zeros(nz); 

%Total opacity 
X_tot = zeros(nz,1); 
%Photon destruction probability eps = X_abs/X_tot
%This is set to 1, as absorption and emission dominate in LTE
destruction_probability = ones(nz,1); 
%Thermal blackbody source function
T = 1000; % in K
temperatures = T*ones(nz,1); %uniform domain temperature
%Isotropic blackbody function
%B_v(T) = 2*h*v^3/c^2 * 1/(e^(hv/kT) -1)
B = 1e-3;
thermal_source = B*temperatures;

%Calculate source function
%Integrate I for mean intensity, J using quadrature rule
%mean_intensity = 1/(4*pi)*sum(weights*ones(1,ntheta)*intensity);
source_function = destruction_probability.*thermal_source; % + ...
%(ones(nz,1) - destruction_probability).*mean_intensity;

%Optical depth intervals before and after pt: depends on opacticy and location
opt_depth = zeros(2,1);
%Interpolation coefficients
interp_coeff = zeros(3,1);

%Test Problems:
%Coherent searchlight beam from LHS of domain
%intensity(1,:) = 1;
%X_tot = zeros(nz,1); %this implies no blackbody emission or absorption

%Uniform material
X_tot = 10*ones(nz,1);

for j=1:ntheta
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
        opt_depth(1) = dz/abs(mu(j))*(X_tot(k)+X_tot(k+upwind)/2);
        %path integral to next point
        opt_depth(2) = dz/abs(mu(j))*(X_tot(k+downwind)+X_tot(k)/2);
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
