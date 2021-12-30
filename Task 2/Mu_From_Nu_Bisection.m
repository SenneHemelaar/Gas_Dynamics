function [mu, mach_i] = Mu_From_Nu_Bisection(nu, mach_left, mach_right, tolerance)
%MU_FROM_NU_BISECTION 
    while (abs(mach_left - mach_right) >= tolerance)
        mach_i = (mach_left + mach_right)/2;
        step = Mach_For_Prandtl_Meyer(mach_i, nu);
        if (abs(step) < abs(tolerance)) 
            break
        elseif (step * Mach_For_Prandtl_Meyer(mach_left, nu) < tolerance)
            mach_right = mach_i;
        else
            mach_left = mach_i;
            mu = Mu_From_Mach(mach_i);
        end
    end
end

