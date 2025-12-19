function [vf, SLt, SLc, STt, STc, SLTs] = getMaterialStrengths()
% getMaterialStrengths - Returns strength properties for Tsai-Wu failure criterion
%
% INPUT:
%   materialChoice - integer (1 to 5), same as used in getMaterialProperties
%
% OUTPUTS:
%   vf   - fiber volume fraction
%   SLt  - Longitudinal tensile strength (ksi)
%   SLc  - Longitudinal compressive strength (ksi)
%   STt  - Transverse tensile strength (ksi)
%   STc  - Transverse compressive strength (ksi)
%   SLTs - In-plane shear strength (ksi)
material_choice = input('Enter the material number again (1-5): ');

switch material_choice
    case 1
        vf   = 0.70;
        SLt  = 217.50;  % ksi
        SLc  = 217.50;
        STt  = 5.80;
        STc  = 35.7;
        SLTs = 9.86;

    case 2
        vf   = 0.50;
        SLt  = 182.70;
        SLc  = 362.50;
        STt  = 8.85;
        STc  = 29.30;
        SLTs = 9.72;

    case 3
        vf   = 0.66;
        SLt  = 209.90;
        SLc  = 209.90;
        STt  = 7.50;
        STc  = 29.90;
        SLTs = 13.50;

    case 4
        vf   = 0.45;
        SLt  = 154;
        SLc  = 88.50;
        STt  = 4.50;
        STc  = 17.10;
        SLTs = 10.40;

    case 5
        vf   = 0.60;
        SLt  = 203;
        SLc  = 34.10;
        STt  = 1.74;
        STc  = 7.69;
        SLTs = 4.93;

    otherwise
        error('Invalid material choice. Please select a number between 1 and 5.');
end

% --- Display selected strengths ---
fprintf('Fiber volume fraction (vf): %.2f\n', vf);
fprintf('Longitudinal tensile strength (SLt): %.2f ksi\n', SLt);
fprintf('Longitudinal compressive strength (SLc): %.2f ksi\n', SLc);
fprintf('Transverse tensile strength (STt): %.2f ksi\n', STt);
fprintf('Transverse compressive strength (STc): %.2f ksi\n', STc);
fprintf('In-plane shear strength (SLTs): %.2f ksi\n\n', SLTs);

end
