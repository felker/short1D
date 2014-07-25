function [interp_coeff] = sc_interpolation(opt_depth_u, opt_depth_d)
%Input:
% opt_depth_u: upwind depth
% opt_depth_d: downwind depth
% Can't have 0 input, leads to NaN via division by zero

interp_coeff = zeros(3,1);

%linear interpolation
%From Kunasz,Auer88
%Is this what they did for their search beam test to avoid NaN?)
if (opt_depth_u ~= 0)
interp_coeff(1) = int_interp_fn(0,opt_depth_u) - int_interp_fn(1,opt_depth_u)/opt_depth_u;
interp_coeff(2) = int_interp_fn(1,opt_depth_u)/opt_depth_u; 
interp_coeff(3) = 0;
end
end

function fn_value = int_interp_fn(i,x)
%Functions used to evaluate the optical depth source function integral
%From Kunasz,Auer88
switch i
    case 0
        fn_value = 1-exp(-x);
    case 1
        fn_value = x -1 + exp(-x);
    case 2
        fn_value = x^2 - 2*x+2*exp(-x);
    otherwise
        error('not a valid interpolating function');
end
end