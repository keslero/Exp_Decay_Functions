% clc, clear, close all
function main_lambdas(Re)

%% Parameters
% Re = 1760;
fle_name = ['lambda/Re_',num2str(Re),'/'];

tm = tm_range(Re);

domains = size(tm,1);
vel_comp = 3; % Velocity components
N = domains*vel_comp;

%% ADD PATHS
% Save location
save_loc = ...
    strcat('/home/shaul/OK/O.A/Exp_Decay/',fle_name);
% Load location
load_loc = ...
    strcat('/media/shaul/Seagate Hub/Omer/Exp_Decay_Data/Re_',num2str(Re),'/');

addpath(['lambda/Re_',num2str(Re)])
for i = 1:length(tm)
    addpath(strcat(load_loc,'POD_SNP-',num2str(tm(i,1)),'-',num2str(tm(i,2)),'/'))
end
addpath('Functions/')
addpath('Figures/')

%% INITIALIZE
fprintf('\nInitializing\n\n')
POD_ID = strcat('POD_SNP-',num2str(tm(1,1)),'-',num2str(tm(1,2)),'/');
filename_r = 'POD200R.cdf.dat';
ncid = netcdf.open([load_loc,POD_ID,filename_r]);        %,'NC_NOWRITE');
[~, n_snap] = netcdf.inqDim(ncid,3);  
netcdf.close(ncid)

%% RETRIEVE APPROXIMATE VELOCITY FIELDS
fprintf('Approximate velocities\n')


lambda_r = zeros(domains,n_snap);
lambda_t = lambda_r;
lambda_z = lambda_r;

for i = 1:domains
    fprintf('\nTime domain %d-%d\n',tm(i,1),tm(i,2))
    [lambda_i] = extract_lam(tm(i,:),1,n_snap,load_loc);

    lambda_r(i,:) = lambda_i(1:n_snap,1);
    lambda_t(i,:) = lambda_i(1:n_snap,2);
    lambda_z(i,:) = lambda_i(1:n_snap,3);
end

%% Save lambdas
% u_r
l = lambda_r(1,:); save(strcat(save_loc,"lambda_r_100_125"),"l");
l = lambda_r(2,:); save(strcat(save_loc,"lambda_r_125_150"),"l");
l = lambda_r(3,:); save(strcat(save_loc,"lambda_r_150_175"),"l");
l = lambda_r(4,:); save(strcat(save_loc,"lambda_r_175_200"),"l");
% u_t
l = lambda_t(1,:); save(strcat(save_loc,"lambda_t_100_125"),"l");
l = lambda_t(2,:); save(strcat(save_loc,"lambda_t_125_150"),"l");
l = lambda_t(3,:); save(strcat(save_loc,"lambda_t_150_175"),"l");
l = lambda_t(4,:); save(strcat(save_loc,"lambda_t_175_200"),"l");
% u_z
l = lambda_z(1,:); save(strcat(save_loc,"lambda_z_100_125"),"l");
l = lambda_z(2,:); save(strcat(save_loc,"lambda_z_125_150"),"l");
l = lambda_z(3,:); save(strcat(save_loc,"lambda_z_150_175"),"l");
l = lambda_z(4,:); save(strcat(save_loc,"lambda_z_175_200"),"l");

%% EIGENSPECTRA
% 
% % Normalizing to first eigenvalue
% lambda_r = lambda_r./lambda_r(:,1);
% lambda_t = lambda_t./lambda_t(:,1);
% lambda_z = lambda_z./lambda_z(:,1);
% 
% % Plot
% % Create [vel_comp X domains] subplots
% % Each column is for a certain snp domain
% % Each row is for a certain vel comp
% 
% % Initiate plot parameters
% fig_ID = figure();
% 
% % Plot
% for i = 1:domains
%     tte = ['SNPs ',num2str(tm(i,1)),' - ',num2str(tm(i,2))];
%     enlarge = 0; % 0 - All eigenvalues | 1 Enlarged area- pick bottom energy limit
%     epsilon = 1e-5; % Lower limit of energy to plot [%]
% 
%     % ~~~~~ Subplot ~~~~~
%     % u_r
%     s = lambda_r(i,:);
%     subplot(vel_comp,domains,1+(i-1))
%     semilogy(s);
%     xlabel('n','FontSize',14)
%     ylabel('\lambda_r','FontSize',16)
%     title(tte);
%     switch enlarge
%         case 1
%             energy = s/sum(s); % Relative energy of each
%             x_max = max(find(energy>epsilon));
%             grid on, xlim([1 x_max]), ylim([epsilon 1])
%         case 0
%             grid on, axis tight
%     end
% 
%     % u_t
%     s = lambda_t(i,:);
%     subplot(vel_comp,domains,1+domains+(i-1))
%     semilogy(s);
%     xlabel('n','FontSize',14)
%     ylabel('\lambda_{\theta}','FontSize',16)
%     title(tte);
%     switch enlarge
%         case 1
%             energy = s/sum(s); % Relative energy of each
%             x_max = max(find(energy>epsilon));
%             grid on, xlim([1 x_max]), ylim([epsilon 1])
%         case 0
%             grid on, axis tight
%     end
% 
%     % u_z
%     s = lambda_z(i,:);
%     subplot(vel_comp,domains,1+2*domains+(i-1))
%     semilogy(s);
%     xlabel('n','FontSize',14)
%     ylabel('\lambda_z','FontSize',16)
%     title(tte);
%     switch enlarge
%         case 1
%             energy = s/sum(s); % Relative energy of each
%             x_max = max(find(energy>epsilon));
%             grid on, xlim([1 x_max]), ylim([epsilon 1])
%         case 0
%             grid on, axis tight
%     end
% 
% end
% 
% % ~~~~~ Save figure ~~~~~
% fig_ID.WindowState = 'maximized';
% switch enlarge
%     case 1
%         saveas(fig_ID,strcat(loc,'Figures/lambdas_enlarged.fig'))
%         saveas(fig_ID,strcat(loc,'Figures/lambdas_enlarged.png'))
%     case 0
%         saveas(fig_ID,strcat(loc,'Figures/lambdas.fig'))
%         saveas(fig_ID,strcat(loc,'Figures/lambdas.png'))
% end

end

%% Functions
function [tm] = tm_range(Re)
switch Re
    case 1760
        tm =  ...
            [100 125
             125 150
             150 175
             175 200];

    case 1780
        tm =  ...
            [125 150
             150 175
             175 200
             200 225];
    case 1800
        tm = [];
    case 1820
        tm =  ...
            [425 450
             450 475
             475 500
             500 525];
end
end
