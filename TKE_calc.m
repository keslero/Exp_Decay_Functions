function [Enrg] = TKE_calc(u_r,u_t,u_z,Nt,Nz,r,choice)
% Calculates energy for all snapshots

%% INITITIALIZE INTEGRATION PARAMETERS
r = flip(1-r)'; % Laminar profile?
r0 = r(1:end-1); r1 = r(2:end);    
R = [0 (r0 + r1)/2 1];
DR = diff(R);
dTh = 2*pi/Nt;
DTh = dTh*ones(1,Nt);  
dA = r.*DR.* DTh';
DA = permute(repmat(dA,1,1,Nz),[2,1,3]); % [r,th,z]
% laminar = (repmat((1-r.^2)',1,Nt,Nz));

%% DEFINE ENERGY 2 be CALCULATED
switch choice
    case 1 % 'r' - Radial
        u_t = 0*u_t;
        u_z = 0*u_z;
    case 2 % 't' - Azimuthal
        u_r = 0*u_r;
        u_z = 0*u_z;
    case 3 % 'z' - Axial
        u_r = 0*u_r;
        u_t = 0*u_t;
    case 4 % 'r-t' - Transverse
        u_z = 0*u_z;
    case 5 % 'r-t-z' - Total
        % Total energy in pipe
end

%% FIND INDICES
Enrg_mat = (u_r.^2+u_t.^2+u_z.^2).*DA; % Transv. energy at all points
Enrg = squeeze(sum(Enrg_mat,[1,2,3])); % Spatial sum of energy

end