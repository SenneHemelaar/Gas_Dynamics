function [mach_for_prandtl_meyer] = Mach_For_Prandtl_Meyer(mach,nu)
%MACH_FOR_PRANDTL_MEYER
mach_for_prandtl_meyer = Prandtl_Meyer_From_Mach(mach) - nu;
end

