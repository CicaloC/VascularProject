function [sigma_extra] = constitutive_model_NH(F)

%=========================================================================
% This function calculates the Cauchy stress extra (sigma_extra) for a  
% given gradient tensor F, Using Neo-Hookean Constitutive model (see Lecture 8)
%=========================================================================

C1 = 80; % material constant, kPa
C = F'*F;

sigma_extra = 2.*C.*C1; % Cauchy extra stress

end


