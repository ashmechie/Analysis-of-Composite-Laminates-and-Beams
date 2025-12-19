function [F1, F2, F11, F22, F66, F12] = getTsaiWuCoefficients(SLt, SLc, STt, STc, SLTs)
% getTsaiWuCoefficients - Calculates Tsai-Wu failure criterion coefficients
%
% INPUTS:
%   SLt  - Longitudinal tensile strength (ksi)
%   SLc  - Longitudinal compressive strength (ksi)
%   STt  - Transverse tensile strength (ksi)
%   STc  - Transverse compressive strength (ksi)
%   SLTs - In-plane shear strength (ksi)
%
% OUTPUTS:
%   F1, F2, F11, F22, F66, F12 - Tsai-Wu coefficients

% --- Tsai-Wu coefficients ---
F11 = 1 / (SLt * SLc);
F22 = 1 / (STt * STc);
F12 = -0.5 * sqrt(F11 * F22);
F66 = 1 / (SLTs^2);
F1  = (1 / SLt) - (1 / SLc);
F2  = (1 / STt) - (1 / STc);

end
