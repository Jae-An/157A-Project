function [MOI] = get_MOI_Cylinder(weight, d_o, d_i, length)
% Gets moment of inertia of a cylindrical shell about its own CG.
% Units: MOI [lb*ft2], weight [lb], length [ft], diameter [ft]

r_o = d_o / 2;
r_i = d_i / 2;
MOI = (1/12) * weight .* (3*(r_o.^2 + r_i.^2) + length.^2);