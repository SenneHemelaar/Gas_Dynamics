function [mach_i, phi_i, nu_i] = Bisection_Method(mach_left, mach_right, tolerance, alpha, backward_phi, backward_nu)
%BISECTION_METHOD Bisection iteration method to find variables (mach, phi
%and nu) in simple wave region ABC
mach_i = (mach_left+mach_right)/2;
while (abs(mach_left - mach_right) >= tolerance)
    mach_i = (mach_left+mach_right)/2;
    step = Bisection_Steps(mach_i, alpha, backward_phi, backward_nu);
    if (abs(step) < abs(tolerance))
        break
    elseif (step * Bisection_Steps(mach_left, alpha, backward_phi, backward_nu) < 0)
        mach_right = mach_i;
    else
        mach_left = mach_i;
    end
end
mu_i  = Mu_From_Mach(mach_i);
phi_i = mu_i + alpha;
nu_i  = Prandtl_Meyer_From_Mach(mach_i);
end