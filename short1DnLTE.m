%1D short characteristics solver for RT
% This is for non-LTE problems using ALI
% Multiple frequencies
%Time independent problem
%Main ref: Hubeny 2003

%Parameters
clear all;
close all;
nz = 1000;
ntheta = 2; %quadarture is only defined up to 12 in each direction, must be even
%order of the quadrature, N, refers to the number of mu-levels in the interval [-1, 1].
nf = 1; %number of frequencies 
MAX_ITER = 10000; %Maximum ALI iterations
delta_c = 1e-4; %Convergence criterion

lz = 10.0;
c = 1.0;
dz = lz/nz;

%Angular Discretization, Sn discrete ordinates. 
[mu, w] = angular_quad1D(ntheta);
mu = [1,-1]';
%Spatial Discretization
%Radiation points are centered on the fluid cells
zz=linspace(0,lz,nz)';

%Monochromatic specific intensity, boundary conditions at 2,nz-1 
intensity = zeros(nz,ntheta,nf); 
%Monochromatic mean intensity
mean_intensity_f = zeros(nz,nf); 
%scattering integral
mean_intensity = zeros(nz,1);

%Total monochromatic opacity 
X_tot = zeros(nz,nf); 
%Photon destruction probability eps = X_abs/X_tot
destruction_prob = ones(nz,1); 
%Thermal blackbody source function
T = 1000; % in K
temperatures = T*ones(nz,1); %uniform domain temperature
%Isotropic blackbody function
%B_v(T) = 2*h*v^3/c^2 * 1/(e^(hv/kT) -1)
B = 1e-3;
thermal_source = B*temperatures;

%Test Problems:
%-------------
%Coherent searchlight beam from LHS of domain
%intensity(1,:,:) = 1;
%X_tot = zeros(nz,nf); %this implies no blackbody emission or absorption

%Uniform material, LTE
%X_tot(:,:) = 1*ones(nz,nf);
%destruction_prob = ones(nz,1);  

%Davis12, highly scattering non-LTE atmosphere
X_tot(:,1) = 50*exp(zz./lz); %exponentially thickening atmosphere

destruction_prob = 1e-2*ones(nz,1);  

%end test problems

%Optical depth intervals using trapezoidal quadrature rule, X(i-1) + X(i)
%Hayek10 computes the optical depth intervals using Bezier interp
%Davis12 simply uses trapezoidal rule
opt_depth = dz*(circshift(X_tot,+1) + X_tot)/2; %1st point is garbage
%Interpolation coefficients
interp_coeff = zeros(3,1);

%Construct the approximate operator
%Ref: Olson, Kunasz 1987
%formal solution, diagonal elements only (we are only taking the i-th index 
%for each i formal solution). Could do tridiagonal, etc
%For linear interpolation, this is a simple formula
normalized_intensity = zeros(nz,ntheta,nf);
for i=1:nf
    for j=1:ntheta/2 %these formulas should be dependent on the sign of mu...
        for k=2:nz-1
            %perhaps because they only depend on the abs value of mu
            %we can ignore the direction, do half of them, assuming
            %symmetry of discretization
            outward_norm = sc_interpolation(opt_depth(k,i)/abs(mu(j)), 0); 
            inward_norm  = sc_interpolation(opt_depth(k+1,i)/abs(mu(j)), 0); 
            normalized_intensity(k,j,i) = 0.5*(outward_norm(2) + inward_norm(2));
        end
       inward_norm  = sc_interpolation(opt_depth(1,i)/abs(mu(j)), 0);
       outward_norm  = sc_interpolation(opt_depth(nz,i)/abs(mu(j)), 0); 
       normalized_intensity(1,j,i) = inward_norm(2); %correct way to handle endpts?
       normalized_intensity(nz,j,i) = outward_norm(2);
    end
end
%Outward going intensity bc, I^+ = I(lz,mu negative)
%Inward going intensity bc, I^- = I(0,mu positive)
%I(0, mu neg... and I(lz, mu positive are free parameters

%integrate over frequnecy and mu
lambda_star = zeros(nz,1);
%frequency weights??
w_f = 1/nf*ones(nf,1);

for i=1:nz
    %should only integrate half the angles?
    lambda_star(i) = (w(1:ntheta/2)'*squeeze(normalized_intensity(i,1:ntheta/2,:)))*w_f;
end

%Initial guess of source function for iterative scheme
source_function = thermal_source; 
for iter=1:MAX_ITER
for i=1:nf 
    for j=1:ntheta
        if mu(j) >=0
            first = 2;
            last = nz;
            upwind = -1; %points to upwind meshpoint
            downwind = 1;
            depth_dir = 0; %points to correct opt_depth interval index
        else
            first = nz-1;
            last = 1;
            depth_dir = 1;
            upwind = 1;
            downwind = -1;
        end
        for k=first:downwind:last %trace characteristics downwind
            %Calculate interpolation coefficients
            %Possible schemes: linear, parabolic, Bezier with control points
            interp_coeff = sc_interpolation(opt_depth(k+depth_dir,i)/abs(mu(j)),opt_depth(k,i)/abs(mu(j)));
            %Update %line above must be fixed for downwind (and below) for
            %parabolic
            intensity(k,j,i) = intensity(k+upwind,j)*exp(-opt_depth(k+depth_dir,i)/abs(mu(j))) + interp_coeff(1)*source_function(k+upwind) + ...
            +interp_coeff(2)*source_function(k) + interp_coeff(3)*source_function(k);
        end
    end
end
%Calculate source function
%Frequency independent for two-level atom (not true... see davis)
%Integrate I for monocrhomatic mean intensity, J_v using quadrature rule
for i=1:nf
    for j=1:nz
        mean_intensity_f(j,i) = 1/2*squeeze(intensity(j,:,i))*w;
    end
end
mean_intensity = sum(mean_intensity_f,2);
 %Update source function
 delta_source = 1./(1 - (1 - destruction_prob).*lambda_star).* ...
     ((1- destruction_prob).*mean_intensity + destruction_prob.*thermal_source - source_function);
 source_function = source_function + delta_source;
 %Convergence check
 iter
 max_rel_change = max(abs(delta_source)./source_function)
 if max_rel_change <= delta_c
     break;
 end
end

%Output

%Analytic solution for exponential scattering atmosphere
for i=2:nz
    optical_depth(i) = sum(opt_depth(2:i,1),1);
end
optical_depth = optical_depth';
analytic_J = thermal_source.*(1-(exp(-sqrt(3).*destruction_prob.*optical_depth))./(1+sqrt(destruction_prob)));
analytic_error = abs(analytic_J - mean_intensity);
loglog(optical_depth,source_function./thermal_source,optical_depth,analytic_J)
%almost matches solution for :
%X_tot(:,1) = 50*exp(zz./lz); %exponentially thickening atmosphere
%destruction_prob = 1e-2*ones(nz,1);  


%plot(zz,analytic_error); 
%h = figure;
%loglog(optical_depth,source_function./thermal_source); %Davis12 Fig 4.
%h2 = figure;
%plot(zz,destruction_prob.*thermal_source,zz,mean_intensity);
%legend('thermal source function','mean intensity');

%Plot each angular intensity (recall, ntheta must be even)
%for i=1:ntheta
%    hi = subplot(2,1,i); 
%    plot(hi,zz,intensity(:,i,1));
%    x_label = sprintf('mu = %f',mu(i));
%    xlabel(x_label);
%end
%h = figure;
%quiver(zeros(ntheta,1),zeros(ntheta,1),mu,sqrt(ones(ntheta,1) - mu.^2)); 
