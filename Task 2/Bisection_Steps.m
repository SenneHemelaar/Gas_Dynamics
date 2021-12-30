function [step] = Bisection_Steps(mach, alpha, backward_phi, backward_nu)
%BISECTION_STEPS
step = - Prandtl_Meyer_From_Mach(mach) + Mu_From_Mach(mach) + alpha - backward_phi + backward_nu;
end

