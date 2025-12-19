function Qbar = calcQbar(Q, theta_deg)
% ===============================================================
% Function to compute transformed reduced stiffness matrices [Qbar]
% for each ply in a laminate.
%
% Inputs:
%   Q         - (3x3) reduced stiffness matrix in material coordinates
%   plyAngles - (1xN) vector of fiber angles (in degrees) for each ply
%
% Outputs:
%   Qbar - (3x3xN) array, where Qbar(:,:,k) is the Qbar matrix of ply k
% ===============================================================

numPlies = length(theta_deg);
Qbar = zeros(3,3,numPlies);   % Preallocate numeric 3D array

for k = 1:numPlies
    theta = theta_deg(k);
    
    % Call existing functions
    T = Tmatrix(theta);
    Tcap = Tcapmatrix(theta);
    % Compute Qbar
    Qbar(:,:,k) = inv(T) * Q * Tcap;
end

end
