function [u_r_out,u_t_out] = puff_finder(u_r,u_t,Nt,Nz,m,r,L,D_range)
% For all snapshots
u_r = squeeze(u_r);
u_t = squeeze(u_t);

% Initialize integration parameters
r = flip(1-r)'; % Laminar profile?
r0 = r(1:end-1); r1 = r(2:end);    
R = [0 (r0 + r1)/2 1];
DR = diff(R);
dTh = 2*pi/Nt;
DTh = dTh*ones(1,Nt);  
dA = r.*DR.* DTh';
DA = permute(repmat(dA,1,1,Nz),[2,1,3]); % [r,th,z]
% laminar = permute (repmat((1-r.^2)',1,Nt,Nz),[3 2 1]); % [z,th,r]

% Initialize matrices
First = 1;
Last  = m;
D = Nz/L; % L is # of diameters in pipe length | Nz is # points in pipe length | Nz/L is # of points per 1 diameter
windLR = round(D_range*D); % Length of puff interval from each side of position
                     % of max TKE
Nzn = (2*windLR+1);  % Total length of puff interval
puff_indices = zeros(m,Nzn);

%% FIND INDICES
TrsEnrg_mat = (u_r.^2+u_t.^2).*DA;   % Trans energy at all points
TrsEnrg = sum(sum(TrsEnrg_mat,2),1); % Trans energy at each cross-sec
TrsEnrg = squeeze(TrsEnrg);

u_r_out = zeros(length(r),Nt,Nzn,m);
u_t_out = zeros(length(r),Nt,Nzn,m);
for i = First:Last
    TrsEnrg_i = squeeze(TrsEnrg(:,i));
    s_max = find(TrsEnrg_i == max(TrsEnrg_i)); % Sec for whitch max TrsEnrg
    sInd = [s_max-windLR : s_max s_max+1 : s_max+windLR];
    % Apply periodicity
    sInd(sInd > Nz | sInd <0) = mod(sInd(sInd > Nz | sInd <0),Nz);
    sInd(sInd == 0) = Nz;

    % puff_indices(i,:) = sInd;

    u_r_i = u_r(:,:,sInd,i);
    u_t_i = u_t(:,:,sInd,i);
    u_r_out(:,:,:,i) = u_r_i;
    u_t_out(:,:,:,i) = u_t_i;
end
end