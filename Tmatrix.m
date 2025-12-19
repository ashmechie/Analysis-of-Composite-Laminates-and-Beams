function T = Tmatrix(theta_deg)
% TMATRIX  Returns the stress transformation matrix [T]
% Input:  theta_deg - ply orientation angle in degrees
% Output: T - 3x3 transformation matrix for stress

theta = deg2rad(theta_deg); % Convert to radians
c = cos(theta);
s = sin(theta);

T = [ c^2      s^2       2*s*c;
      s^2      c^2      -2*s*c;
     -s*c      s*c       (c^2 - s^2) ];
end

