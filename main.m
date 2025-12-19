%% Final main script: progressive failure (zero stiffness), store R at failure per ply,
% then pick the maximum of those stored R's to compute critical loads.

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
Exx = (A(1,1)/h_total) * (1 - (A(1,2)^2)/(A(1,1)*A(2,2)));
fprintf('\nEffective Exx of laminate = %.6e psi\n', Exx);

% ---- NEW: Compute effective laminate Gxy ----
Gxy_eff = A(3,3) / h_total;
fprintf('Effective Gxy of laminate = %.6e psi\n', Gxy_eff);
if rcond(ABD) < 1e-14
    error('Full-laminate ABD is singular or ill-conditioned.');
end

strain_curvature = ABD \ loads;
eps0 = strain_curvature(1:3);
kappa = strain_curvature(4:6);

%% ---------- Compute and Display Ply Stresses/Strains (Local) ----------
eps_local_all = zeros(3,n);
sig_local_all = zeros(3,n);

for i = 1:n
    strain_global = eps0 + z_mid(i) * kappa;
    Tcap = Tcapmatrix(theta_deg(i));
    eps_local_all(:,i) = Tcap * strain_global;
    sig_local_all(:,i) = Q * eps_local_all(:,i);
end


%% ---------- Initial R-values ----------
R_initial = NaN(n,1);
for i = 1:n
    s11 = sig_local_all(1,i)/1000;
    s22 = sig_local_all(2,i)/1000;
    t12 = sig_local_all(3,i)/1000;

    Acoef = F11*s11^2 + 2*F12*s11*s22 + F22*s22^2 + F66*t12^2;
    Bcoef = F1*s11 + F2*s22;
    Ccoef = -1;
    disc = Bcoef^2 - 4*Acoef*Ccoef;
    if disc >= 0
        R1 = (-Bcoef + sqrt(disc)) / (2*Acoef);
        R2 = (-Bcoef - sqrt(disc)) / (2*Acoef);
        R_initial(i) = max([R1,R2,0]);
    else
        R_initial(i) = NaN;
    end
end

fprintf('\nInitial R-values (full laminate):\n');
disp(table((1:n)', theta_deg', z_mid, R_initial, ...
    'VariableNames', {'Ply_Number','Angle_deg','z_mid_in','R_initial'}));

%% ---------- Progressive Failure Loop (NO per-iteration tables) ----------
active = true(n,1);
R_fail_at_ply = NaN(n,1);
failure_order = [];     % track order of failures
iteration = 1;
tol = 1e-9;

fprintf('\nStarting progressive failure loop...\n');

while true
    % Apply zero stiffness to failed plies
    Qbar_mod = zeros(3,3,n);
    for i = 1:n
        if active(i)
            Qbar_mod(:,:,i) = Qbar_all(:,:,i);
        else
            Qbar_mod(:,:,i) = zeros(3);
        end
    end
    
    % Recompute ABD
    A = zeros(3); B = zeros(3); D = zeros(3);
    for i = 1:n
        z1 = z_interfaces(i); z2 = z_interfaces(i+1);
        Qb = Qbar_mod(:,:,i);
        A = A + Qb * (z2 - z1);
        B = B + 0.5 * Qb * (z2^2 - z1^2);
        D = D + (1/3) * Qb * (z2^3 - z1^3);
    end
    ABD = [A B; B D];
    
    if rcond(ABD) < 1e-12
        warning('ABD singular or nearly singular. Stopping progressive failure loop.');
        break;
    end
    
    strain_curvature = ABD \ loads;
    eps0 = strain_curvature(1:3);
    kappa = strain_curvature(4:6);
    
    % --- Compute new R-values
    R_iter = NaN(n,1);
    for i = 1:n
        strain_global = eps0 + z_mid(i) * kappa;
        Tcap = Tcapmatrix(theta_deg(i));
        strain_local = Tcap * strain_global;
        stress_local = Q * strain_local;
        
        s11 = stress_local(1)/1000;
        s22 = stress_local(2)/1000;
        t12 = stress_local(3)/1000;
        Acoef = F11*s11^2 + 2*F12*s11*s22 + F22*s22^2 + F66*t12^2;
        Bcoef = F1*s11 + F2*s22;
        Ccoef = -1;
        disc = Bcoef^2 - 4*Acoef*Ccoef;
        if disc >= 0
            R1 = (-Bcoef + sqrt(disc)) / (2*Acoef);
            R2 = (-Bcoef - sqrt(disc)) / (2*Acoef);
            R_iter(i) = max([R1,R2,0]);
        end
    end
    
    % --- Identify ply to fail
    active_indices = find(active);
    active_Rs = R_iter(active_indices);
    valid_mask = ~isnan(active_Rs);
    if ~any(valid_mask)
        fprintf('No valid R among active plies. Stopping.\n');
        break;
    end
    
    minR = min(active_Rs(valid_mask));
    to_fail_mask = (R_iter <= minR + tol) & active;
    to_fail_indices = find(to_fail_mask);
    if isempty(to_fail_indices)
        fprintf('No plies identified to fail this iteration. Stopping.\n');
        break;
    end
    
    % --- Store failure info
    for idx = to_fail_indices'
        R_fail_at_ply(idx) = minR;
        failure_order = [failure_order; iteration, idx, minR];
    end
    
    % Mark ply(s) as failed
    active(to_fail_indices) = false;
    
    % End condition
    if ~any(active)
        fprintf('All plies have failed. End of progressive failure.\n');
        break;
    end
    
    iteration = iteration + 1;
end

%% ---------- Ply Failure Order Summary (NEW TABLE) ----------
if ~isempty(failure_order)
    fprintf('\n=============================================\n');
    fprintf('Ply Failure Order Summary:\n');
    FailureSummary = array2table(failure_order, ...
        'VariableNames', {'Iteration','Ply_Number','R_fail'});
    disp(FailureSummary);
else
    fprintf('\nNo ply failures recorded.\n');
end
fprintf('=============================================\n');


%% ---------- Post-Failure Analysis ----------
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
%% ---------- Compute and Display Ply Stresses/Strains (Local & Global) ----------
eps_local_all = zeros(3,n);
sig_local_all = zeros(3,n);
eps_global_all = zeros(3,n);
sig_global_all = zeros(3,n);

fprintf('\nPly-wise Local Strains and Stresses (before R-value calculation):\n');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
fprintf('Ply |  Angle(deg) |        ε1            ε2            γ12      |        σ1 (psi)        σ2 (psi)        τ12 (psi)\n');
fprintf('-------------------------------------------------------------------------------------------------------------\n');

for i = 1:n
    % ---- Global strain (εx, εy, γxy) ----
    strain_global = eps0 + z_mid(i) * kappa;
    eps_global_all(:,i) = strain_global;

    % ---- Transform to local coordinates ----
    Tcap = Tcapmatrix(theta_deg(i));
    eps_local_all(:,i) = Tcap * strain_global;
    sig_local_all(:,i) = Q * eps_local_all(:,i);

    % ---- Compute global stresses (σx, σy, τxy) ----
    % Note: use transformation from local to global
    Tinv = inv(Tcap);   % since [ε_local] = [Tcap]*[ε_global]
    sig_global_all(:,i) = Tinv' * sig_local_all(:,i);  

    % ---- Print local values ----
    fprintf('%3d | %10.2f | %12.6e %12.6e %12.6e | %14.3f %14.3f %14.3f\n', ...
        i, theta_deg(i), ...
        eps_local_all(1,i), eps_local_all(2,i), eps_local_all(3,i), ...
        sig_local_all(1,i), sig_local_all(2,i), sig_local_all(3,i));
end
fprintf('-------------------------------------------------------------------------------------------------------------\n');

%% ---------- NEW: Display Global Strains and Stresses Table ----------
fprintf('\nPly-wise Global Strains and Stresses (x-y coordinates):\n');
fprintf('-------------------------------------------------------------------------------------------------------------\n');
fprintf('Ply |  Angle(deg) |        εx            εy            γxy      |        σx (psi)        σy (psi)        τxy (psi)\n');
fprintf('-------------------------------------------------------------------------------------------------------------\n');

for i = 1:n
    fprintf('%3d | %10.2f | %12.6e %12.6e %12.6e | %14.3f %14.3f %14.3f\n', ...
        i, theta_deg(i), ...
        eps_global_all(1,i), eps_global_all(2,i), eps_global_all(3,i), ...
        sig_global_all(1,i), sig_global_all(2,i), sig_global_all(3,i));
end
fprintf('-------------------------------------------------------------------------------------------------------------\n');

%% ---------- PLOTTING STRESSES AND STRAINS (SUBPLOTS) ----------
z_plot = [];
strains_plot = []; 
stresses_plot = []; 

for i = 1:n
    z_top = z_interfaces(i);
    z_bot = z_interfaces(i+1);
    
    % Using initial Qbar_all for these plots
    Qb = Qbar_all(:,:,i); 
    
    z_points = [z_top, z_bot];
    for zp = z_points
        eps_g = eps0 + zp * kappa;
        sig_g = Qb * eps_g;
        
        z_plot = [z_plot; zp];
        strains_plot = [strains_plot; eps_g'];
        stresses_plot = [stresses_plot; sig_g'];
    end
end

% --- Figure 1: Global Strains Subplots ---
figure('Color', 'w', 'Name', 'Global Strains Breakdown', 'Position', [100, 100, 1000, 400]);

subplot(1,3,1);
plot(strains_plot(:,1), z_plot, 'r', 'LineWidth', 2); grid on;
title('\epsilon_x (Normal X)'); ylabel('z (in)'); xlabel('Strain');

subplot(1,3,2);
plot(strains_plot(:,2), z_plot, 'b', 'LineWidth', 2); grid on;
title('\epsilon_y (Normal Y)'); xlabel('Strain');

subplot(1,3,3);
plot(strains_plot(:,3), z_plot, 'g', 'LineWidth', 2); grid on;
title('\gamma_{xy} (Shear)'); xlabel('Strain');

sgtitle('Global Strains Through Thickness');

% --- Figure 2: Global Stresses Subplots ---
figure('Color', 'w', 'Name', 'Global Stresses Breakdown', 'Position', [100, 550, 1000, 400]);

subplot(1,3,1);
plot(stresses_plot(:,1), z_plot, 'r', 'LineWidth', 2); grid on;
title('\sigma_x (Normal X)'); ylabel('z (in)'); xlabel('Stress (psi)');

subplot(1,3,2);
plot(stresses_plot(:,2), z_plot, 'b', 'LineWidth', 2); grid on;
title('\sigma_y (Normal Y)'); xlabel('Stress (psi)');

subplot(1,3,3);
plot(stresses_plot(:,3), z_plot, 'g', 'LineWidth', 2); grid on;
title('\tau_{xy} (Shear)'); xlabel('Stress (psi)');

sgtitle('Global Stresses Through Thickness');
