function [prandtl_meyer] = Prandtl_Meyer_From_Mach(mach)
%Prandtl_Meyer_From_Mach: Calculates the Prandtl-Meyer in radians from given mach number
gamma = 7/5;
prandtl_meyer = sqrt((gamma+1)/(gamma-1)) * atan( sqrt(((gamma-1)/(gamma+1))...
                * (mach^2-1)) ) - atan(sqrt(mach^2 -1));
end

