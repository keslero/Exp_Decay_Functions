function [lambda] = extract_lam(tm,np,nq,load_loc)
% clc, clear, tm = [25 50]; np = 1; nq = 2;

%% PATHS & FILE NAMES
path = strcat(load_loc,'POD_SNP-', ...
    num2str(tm(1)),'-',num2str(tm(2)),'/');
filename_r = strcat(path,'POD200R.cdf.dat');
filename_t = strcat(path,'POD200T.cdf.dat');
filename_z = strcat(path,'POD200Z.cdf.dat');
stateName = strcat(path,'state',num2str(6112),'.cdf.dat');

%% INITIALIZE
% fprintf('Initializing\n\n')

ncid = netcdf.open(filename_r);      %,'NC_NOWRITE');
[~, Nr]     = netcdf.inqDim(ncid,0); % Radial points
[~, Nt]     = netcdf.inqDim(ncid,1); % Azimutal points
[~, Nz]     = netcdf.inqDim(ncid,2); % Axial points
[~, n_snap] = netcdf.inqDim(ncid,3);  
netcdf.close(ncid)

vel_r = ncread(filename_r,'matVar'); % vel field r [Nz,Nt,Nr,n_snp]
vel_t = ncread(filename_t,'matVar'); % vel field t [Nz,Nt,Nr,n_snp]
vel_z = ncread(filename_z,'matVar'); % vel field z [Nz,Nt,Nr,n_snp]

vel_r = permute(vel_r,[3,2,1,4]);    % vel field r [Nr,Nt,Nz,n_snp]
vel_t = permute(vel_t,[3,2,1,4]);    % vel field t [Nr,Nt,Nz,n_snp]
vel_z = permute(vel_z,[3,2,1,4]);    % vel field z [Nr,Nt,Nz,n_snp]

vel_r = reshape(vel_r,[],n_snap);
vel_t = reshape(vel_t,[],n_snap);
vel_z = reshape(vel_z,[],n_snap);

%% POD
switch 1
    case 1
        % Omer's POD
        fprintf('\n\tOmer''s method\n')

        fprintf('u_r\n')
        [lambda_r,a_r,phi_r] = POD(n_snap,Nr,Nt,Nz,vel_r);
        fprintf('u_t\n')
        [lambda_t,a_t,phi_t] = POD(n_snap,Nr,Nt,Nz,vel_t);
        fprintf('u_z\n')
        [lambda_z,a_z,phi_z] = POD(n_snap,Nr,Nt,Nz,vel_z);

        a_r = a_r(:,np:nq);
        a_t = a_t(:,np:nq);
        a_z = a_z(:,np:nq);

        phi_r = phi_r(np:nq,:);
        phi_t = phi_t(np:nq,:);
        phi_z = phi_z(np:nq,:);

        u_r = a_r*phi_r;
        u_t = a_t*phi_t;
        u_z = a_z*phi_z;

        u_r = reshape(u_r',Nr,Nt,Nz,n_snap);
        u_t = reshape(u_t',Nr,Nt,Nz,n_snap);
        u_z = reshape(u_z',Nr,Nt,Nz,n_snap);

    case 2
        % Alex's POD
        fprintf('\tAlex''s method\n')

        cenergy = 0; % Not relevant to me
        [phi_r,at_r,lambda_r,~,nbasis_r,~] = POD2(vel_r,cenergy);
        [phi_t,at_t,lambda_t,~,nbasis_t,~] = POD2(vel_t,cenergy);
        [phi_z,at_z,lambda_z,~,nbasis_z,~] = POD2(vel_z,cenergy);

        % VELOCITY APPROXIMATION
        fprintf('Truncating velocity fields\n')

        u_r = PODapprox2(at_r,phi_r,n_snap,np,nq);
        u_t = PODapprox2(at_t,phi_t,n_snap,np,nq);
        u_z = PODapprox2(at_z,phi_z,n_snap,np,nq);

        u_r = reshape(u_r,Nr,Nt,Nz,n_snap);
        u_t = reshape(u_t,Nr,Nt,Nz,n_snap);
        u_z = reshape(u_z,Nr,Nt,Nz,n_snap);
end

lambda = [lambda_r lambda_t lambda_z];

end