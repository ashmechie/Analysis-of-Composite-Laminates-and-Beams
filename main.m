%% ------------------------------------------------------------------
% Progressive failure analysis (zero stiffness) for laminates
% Only prints final summary table of R_fail per ply
% ------------------------------------------------------------------

clc; clear; close all;

%% ---------- Material & strengths ----------
[E11, E22, G12, v12, v21] = getMaterialProperties();
[vf, SLt, SLc, STt, STc, SLTs] = getMaterialStrengths();

% Reduced stiffness Q (plane stress)
Q11 = E11 / (1 - v12 * v21);
Q12 = v12 * E22 / (1 - v12 * v21);
Q22 = E22 / (1 - v12 * v21);
Q66 = G12;
Q = [ Q11 Q12 0; Q12 Q22 0; 0 0 Q66 ];

%% ---------- Laminate geometry & stacking ----------
t = input('Enter thickness of each ply (in inches): ');

fprintf('Enter ply orientations as [angle number_of_plies].\n');
fprintf('Example: 45 10 (means 10 plies at 45 deg)\n');
fprintf('Enter 000 to finish entry.\n');

theta_deg = [];  % initialize empty
while true
    ply_input = input('Enter angle and number of plies: ','s');  % read as string
    if strcmp(ply_input, '000')
        break;  % end entry
    end
    
    % Parse the input into two numbers
    nums = sscanf(ply_input, '%f %d');  % angle (float), count (integer)
    if numel(nums) ~= 2
        fprintf('Invalid input. Enter as: angle number_of_plies, e.g., 45 10\n');
        continue;
    end
    
    angle = nums(1);
    count = nums(2);
    
    % Append repeated angles
    theta_deg = [theta_deg, repmat(angle, 1, count)];
end

% Total number of plies
n = length(theta_deg);

% Compute laminate mid-plane z-coordinates
h_total = n * t;
z_interfaces = linspace(-h_total/2, h_total/2, n+1);
z_mid = zeros(n,1);
for i = 1:n
    z_mid(i) = 0.5 * (z_interfaces(i) + z_interfaces(i+1));
end

fprintf('\nTotal plies = %d\n', n);
fprintf('Stacking sequence (top->bottom):\n');
disp(theta_deg);
fprintf('\nPly mid-surface z (in), top->bottom:\n');
disp(z_mid);

%% ---------- Loads from user ----------
% loads = [Nx; Ny; Nxy; Mx; My; Mxy]
loads = inputLoads();
if numel(loads) ~= 6
    error('inputLoads must return a 6-element column vector [Nx;Ny;Nxy;Mx;My;Mxy].');
end
loads = loads(:); % ensure column
load_names = {'Nx','Ny','Nxy','Mx','My','Mxy'};

%% ---------- Tsai-Wu coeffs ----------
[F1,F2,F11,F22,F66,F12] = getTsaiWuCoefficients(SLt, SLc, STt, STc, SLTs);

%% ---------- Initial ABD ----------
Qbar_all = calcQbar(Q, theta_deg);  % expected 3x3xn
if ndims(Qbar_all) == 2
    Qbar_stack = zeros(3,3,n);
    for ii = 1:n
        Qbar_stack(:,:,ii) = calcQbar(Q, theta_deg(ii));
    end
    Qbar_all = Qbar_stack;
end

[A, B, D, ABD] = calcABD(Qbar_all, theta_deg, t);

if rcond(ABD) < 1e-14
    error('Full-laminate ABD is singular or ill-conditioned.');
end

strain_curvature = ABD \ loads;
eps0 = strain_curvature(1:3);
kappa = strain_curvature(4:6);

%% ---------- Call the failure-computation function ----------
[R_initial, R_fail_at_ply, failure_order] = computeInitialAndProgressiveFailure( ...
        n, theta_deg, z_mid, z_interfaces, Q, Qbar_all, ...
        F1, F2, F11, F22, F66, F12, loads, eps0, kappa);

%% ---------- Post-Failure Analysis (summary only) ----------
if all(isnan(R_fail_at_ply))
    error('No ply failure records were captured (R_fail_at_ply all NaN).');
end

R_critical = max(R_fail_at_ply(~isnan(R_fail_at_ply)));
ply_indices_with_Rcrit = find(abs(R_fail_at_ply - R_critical) < 1e-12);
ply_of_Rcrit = ply_indices_with_Rcrit(1);

fprintf('\n=============================================\n');
fprintf('Stored R_fail per ply (NaN = did not record fail):\n');
disp(table((1:n)', theta_deg', z_mid, R_fail_at_ply, ...
    'VariableNames', {'Ply_Number','Angle_deg','z_mid_in','R_fail_recorded'}));

fprintf('\nSelected R_critical (max among stored per-ply failures): %g\n', R_critical);
fprintf('Corresponding ply (first match) = Ply #%d, angle = %g deg, z_mid = %g in\n', ...
    ply_of_Rcrit, theta_deg(ply_of_Rcrit), z_mid(ply_of_Rcrit));

%% ---------- Critical Loads ----------
critical_loads = R_critical * loads;
fprintf('\nCritical loads (R_critical × applied load components):\n');
for k = 1:6
    if loads(k) ~= 0
        fprintf(' %3s : Applied = %g  -->  Critical = %g\n', load_names{k}, loads(k), critical_loads(k));
    end
end
fprintf('=============================================\n');
%% ---------- Stress and Strain Plotting ----------
% Calculate stresses and strains for each ply based on INITIAL state
% (Before progressive degradation)
z_plot = linspace(-h_total/2, h_total/2, 100);
stress_plot = zeros(length(z_plot), 3);
strain_plot = zeros(length(z_plot), 3);

for i = 1:length(z_plot)
    z = z_plot(i);
    % Find which ply this z belongs to
    ply_idx = find(z <= z_interfaces(2:end), 1, 'first');
    if isempty(ply_idx), ply_idx = n; end
    
    % Global Strain at height z: eps = eps0 + z * kappa
    eps_z = eps0 + z * kappa;
    strain_plot(i, :) = eps_z';
    
    % Global Stress at height z: sigma = Qbar * eps
    % Note: Using the original Qbar_all for initial stress state
    sigma_z = Qbar_all(:,:,ply_idx) * eps_z;
    stress_plot(i, :) = sigma_z';
end

% Create Figure for Strains
figure('Name', 'Laminate Strain Distribution', 'NumberTitle', 'off');
labels = {'\epsilon_{xx}', '\epsilon_{yy}', '\gamma_{xy}'};
for i = 1:3
    subplot(3,1,i)
    plot(strain_plot(:,i), z_plot, 'b', 'LineWidth', 1.5)
    grid on; ylabel('z (in)'); xlabel(labels{i});
    title(['Global Strain ', labels{i}]);
end

% Create Figure for Stresses
figure('Name', 'Laminate Stress Distribution', 'NumberTitle', 'off');
labels_sig = {'\sigma_{xx}', '\sigma_{yy}', '\tau_{xy}'};
for i = 1:3
    subplot(3,1,i)
    plot(stress_plot(:,i), z_plot, 'r', 'LineWidth', 1.5)
    grid on; ylabel('z (in)'); xlabel(labels_sig{i});
    title(['Global Stress ', labels_sig{i}]);
end


