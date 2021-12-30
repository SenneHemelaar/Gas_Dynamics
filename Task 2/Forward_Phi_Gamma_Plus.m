function [phi_forward] = Forward_Phi_Gamma_Plus(nu_backward,nu_forward,phi_backward)
%FORWARD_PHI_GAMMA_PLUS Calculates forward phi
phi_forward = nu_forward + phi_backward - nu_backward;
end

