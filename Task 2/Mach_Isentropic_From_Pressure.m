function [mach_2] = Mach_Isentropic_From_Pressure(mach_1,p_1,p_2,gamma)
%MACH_ISENTROPIC_FROM_PRESSURE Computing the mach number using isentropic
%relation
mach_2 = sqrt((p_1/p_2) * (2/(gamma-1)) * (1+ (((gamma-1)/2)*(mach_1^2)) - 1));
end

