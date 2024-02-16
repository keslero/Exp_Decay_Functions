function main(tm,np,nq,choice)
% Inputs:
% -> tm
% -> np
% -> nq
% -> choice (which energies to plot):
%        case 1 % 'r'     - Radial
%        case 2 % 't'     - Azimuthal
%        case 3 % 'z'     - Axial
%        case 4 % 'r-t'   - Transverse
%        case 5 % 'r-t-z' - Total

% Manual operation
% clear, clc
% tm = [100 125
%     175 200
%     200 225];
% np = [1 1 1];
% nq = [1 2 5];
% choice = 1:5;


%% Parameters
domains = size(tm,1);
mode_combos = length(np);
N = domains*mode_combos;

switch mode_combos
    case 1
        figName = sprintf('Modes_%02d-%02d', ...
            np(1),nq(1));
    case 2
        figName = sprintf('Modes_%02d-%02d_%02d-%02d', ...
            np(1),nq(1),np(2),nq(2));
    case 3
        figName = sprintf('Modes_%02d-%02d_%02d-%02d_%02d-%02d', ...
             np(1),nq(1),np(2),nq(2),np(3),nq(3));
end
folder = strcat("Figures/",figName);
mkdir(folder)

%% INITIALIZE
filename_r = 'POD200R.cdf.dat';
filename_t = 'POD200T.cdf.dat';
filename_z = 'POD200Z.cdf.dat';
POD_ID = strcat('POD_SNP-',num2str(tm(1,1)),'-',num2str(tm(1,2)));
loc =  strcat('/home/shaul/OK/O.A/Exp_Decay/');
path =strcat(loc,POD_ID,'/');
stateName = strcat(path,'state6112.cdf.dat');
fprintf('\nInitializing\n\n')
ncid = netcdf.open(filename_r);        %,'NC_NOWRITE');
[~, Nr]     = netcdf.inqDim(ncid,0); % radial points
[~, Nt]     = netcdf.inqDim(ncid,1); % Azimutal points
[~, Nz]     = netcdf.inqDim(ncid,2); % axial points
[~, n_snap] = netcdf.inqDim(ncid,3);  
netcdf.close(ncid)

r = ncread(stateName,'r');
stid = ncinfo(stateName);
para = [stid.Attributes.Value];
alpha = para(3);
L = pi/alpha;

for i = 1:mode_combos
    if nq(i) > n_snap
        nq(i) = n_snap;
    end
end

%% RETRIEVE APPROXIMATE VELOCITY FIELDS
fprintf('Approximate velocities\n')

u_r_app = zeros(N,Nr,Nt,Nz,n_snap); u_t_app = u_r_app; u_z_app = u_r_app;
lambda  = zeros(N,n_snap,3);
k = 1; % Initialize counter
for j = 1:mode_combos
    for i = 1:domains
        fprintf('\n(%d,%d)\n',i,j)
        [u_r_k,u_t_k,u_z_k,lambda_k] = trunc_vels(tm(i,:),np(j),nq(j));

        u_r_app(k,:,:,:,:) = u_r_k(:,:,:,1:n_snap);
        u_t_app(k,:,:,:,:) = u_t_k(:,:,:,1:n_snap);
        u_z_app(k,:,:,:,:) = u_z_k(:,:,:,1:n_snap);
        lambda (k,:,:)     = lambda_k(1:n_snap,:);

        k = k + 1;
    end
end

save(strcat(folder,"/lambda.mat"),"lambda")

%% PUFF FoR
% fprintf('\nChange to puff frame of reference\n')
% 
% clear u_r_k u_t_k
% D_range = 4; % All snp will now be centralized around the max TKE with a domain of +-4D
% for k = 1:N
%     u_r_k = u_app_r(k,:,:,:,:);
%     u_t_k = u_app_t(k,:,:,:,:);
% 
%     [u_r_puff,u_t_puff] = puff_finder(u_r_k,u_t_k,Nt,Nz,n_snap,r,L,D_range);
% 
%     if k == 1 
%             Nz_puff = size(u_r_puff,3);
%             u_r = zeros(N,Nr,Nt,Nz_puff,n_snap);
%             u_t = zeros(N,Nr,Nt,Nz_puff,n_snap);
%     end
% 
%     u_r(k,:,:,:,:) = u_r_puff;
%     u_t(k,:,:,:,:) = u_t_puff;
% end

u_r = u_r_app;
u_t = u_t_app;
u_z = u_z_app;

save(strcat(folder,"/velocities.mat"),"u_r","u_t","u_z")

%% LOOP 4 ALL ENERGY CALCULATIONS
for j = choice

    % ENERGY CALCULATION
    fprintf('\nCalculating energy for all snapshots\n\n')
    TKE = zeros(N,n_snap);
    for k = 1:N
        u_r_k = squeeze(u_r(k,:,:,:,:));
        u_t_k = squeeze(u_t(k,:,:,:,:));
        u_z_k = squeeze(u_z(k,:,:,:,:));
        TKE(k,:) = TKE_calc(u_r_k,u_t_k,u_z_k,Nt,Nz,r,j);
    end

    % PLOT
    % Create [mode_combos X domains] subplots
    % Each column is for a certain snp domain
    % Each row is for a certain mode combination

    switch j
        case 1 % 'r' - Radial
            Enrg_name = 'E_r';
            y_lab     = 'E_r';
        case 2 % 't' - Azimuthal
            Enrg_name = 'E_t';
            y_lab     = 'E_t';
        case 3 % 'z' - Axial
            Enrg_name = 'E_z';
            y_lab     = 'E_z';
        case 4 % 'r-t' - Transverse
            Enrg_name = 'E_trans';
            y_lab = 'E_{\perp}';
        case 5 % 'r-t-z' - Total
            Enrg_name = 'E_tot';
            y_lab = 'E_{tot}';
    end

    % Initiate plot parameters
    tau_D = zeros(N,1); a = tau_D; R2 = tau_D;
    R2_lim = 0.7;
    fig_ID = figure();

    % Plot
    for i = 1:N
        % Find indices
        i_POD  = mod(i,domains); if (i_POD == 0), i_POD = domains; end
        i_mode = floor((i-1)/domains)+1;
        % Create titles as strings
        POD_ID  = strcat(['POD SNP ',num2str(tm(i_POD,1)),'-',num2str(tm(i_POD,2))]);
        mode_ID = strcat(num2str(np(i_mode)),'-',num2str(nq(i_mode)));
        tte = {POD_ID;['Modes ',mode_ID]};

        % Create time vector for current interval
        t = linspace(tm(i_POD,1),tm(i_POD,2),n_snap);
        TKE_i = squeeze(TKE(i,:));

        % Subplot
        subplot(mode_combos,domains,i)
        semilogy(t,TKE_i);

        % Find slope of exp. decay
        go = 1;
        t_slope   = t; 
        TKE_slope = TKE_i;
        while go == 1
            [tau_D_i,a_i] = find_slope(t_slope,TKE_slope);
            f = 10.^(a_i-t_slope/tau_D_i);
            R2_i = coef_det(TKE_slope,f);
            if ( R2_i < R2_lim ) && ( length(t_slope) > n_snap/5 )
                t_slope   = t_slope  (n_snap/5+1:end);
                TKE_slope = TKE_slope(n_snap/5+1:end);
            else
                go = 0;
            end
        end
    
        % If straight enough line (== exp. decay) -> print slope
        if R2_i > R2_lim
            if length(t_slope) == n_snap
                str = sprintf('\\tau_D = %.2f\nR^2 = %.2f',tau_D_i,R2_i);
            else
                str = sprintf( ...
                    '\\tau_D = %.2f\nR^2 = %.2f\n''t = %.0f -> %.0f' ...
                    ,tau_D_i,R2_i,t_slope(1),t_slope(end));
            end
            y_i = 55;
            text(t(y_i+40),TKE_i(y_i),str,"FontSize",12)
        end
        tau_D(i) = tau_D_i;
        a(i) = a_i;
        R2(i) = R2_i;

        % Aesthetics
        grid on, axis tight
        xlabel('t [D/U_{b}]')
        ylabel(strcat(y_lab,' (Domain Total)'))
        title(tte);
    end

    % ~~~~~ Save figure ~~~~~
    fig_ID.WindowState = 'maximized';
    saveas(fig_ID,strcat(loc,folder,'/',Enrg_name,'.fig'))
    saveas(fig_ID,strcat(loc,folder,'/',Enrg_name,'.png'))
    close(fig_ID )

end