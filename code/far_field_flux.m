function [F_far_field] = far_field_flux(W_i, nu_i_inf)

global W_L;
global W_R;

if nu_i_inf(1) < 0
    F_far_field = Roe_2D_flux(W_i, W_L, nu_i_inf);
elseif nu_i_inf(1) > 0
    F_far_field = Roe_2D_flux(W_i, W_R, nu_i_inf);
end


end

