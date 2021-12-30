function [mu] = Mu_From_Mach(mach)
%MU_FROM_MACH Calculates angle mu in radians from mach number
mu = asin(1/mach);
end

