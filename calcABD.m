function [A, B, D, ABD] = calcABD(Qbar, theta_deg, t)
% ===============================================================
% Function to compute the [A], [B], [D] and full [ABD] matrices
% for a composite laminate.
%
% Inputs:
%   Qbar      - (3x3xN) array of Qbar matrices for each ply
%   plyAngles - (1xN) vector of fiber angles
%   t         - thickness of each ply (same for all plies)
%
% Outputs:
%   A   - (3x3) extensional stiffness matrix
%   B   - (3x3) coupling stiffness matrix
%   D   - (3x3) bending stiffness matrix
%   ABD - (6x6) assembled laminate stiffness matrix
% ===============================================================

numPlies = length(theta_deg);

% Laminate total thickness
h_total = numPlies * t;

% Ply interface coordinates through the thickness
z = linspace(-h_total/2, h_total/2, numPlies+1);

% Initialize matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

% Summation over plies
for k = 1:numPlies
    Qk = Qbar(:,:,k);
    z_top = z(k+1);
    z_bot = z(k);
    
    A = A + Qk * (z_top - z_bot);
    B = B + 0.5 * Qk * (z_top^2 - z_bot^2);
    D = D + (1/3) * Qk * (z_top^3 - z_bot^3);
end

% Assemble ABD matrix
ABD = [A B; B D];

end
