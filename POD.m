function [lambda,a,phi] = POD(m,Nr,Nt,Nz,u)
%% Input: 
%         m  = number of snapshots 
%         n  = number of modes (=m)
%         Nr = x-number of points
%         Nt = y-number of points
%         Nz = z-number of points
%         u_r = u_r(1:Nr,1:Nt,1:Nz,m) ur-velocity
%         u_t = u_t(1:Nr,1:Nt,1:Nz,m) ut-velocity
%         u_z = u_z(1:Nr,1:Nt,1:Nz,m) uz-velocity

%% Output: 
%         D = Modal energy
%         FourCoef = temporal coefficients
%         psi,phi,xi = ur-vel POD & ut-vel POD modes & uz-vel POD modes

%% Initialize
n = m;

%% Compute covariance matrix
C = u'*u/m;

%% Compute eigenvalues and eigenvectors of C
[V,lambda] = eig(C);
[lambda,I] = sort(diag(lambda),'descend'); 
V = V(:,I);

%% Scale Fourier Coefficients: <a_i a_i> = lambda_i
% a = zeros(m,n); % These are the Fourier coefficients
% for i = 1:n
%     fprintf('\n\tSingular value %d\n',i)
%     for j = 1:m
%         % ModeFactor = sqrt(m*lambda(i));
%         % a(j,i) = V(j,i)*ModeFactor;
%     end
% end

fprintf('Calculating singular values\n')
ModeFactor = sqrt(m*lambda);
a = V.*ModeFactor';

%% Compute POD modes 
% phi = zeros(size(u,1),n);
% for i = 1:n % Run over modes
%     fprintf('\n\tPOD mode %d\n',i)
%     for j = 1:m % Run over snp
%         phi(:,i) = phi(:,i) + V(j,i)*u(:,i);
%     end
%     % Normalize
%     modeFactor = 1 ./ sqrt(m*lambda(i));
%     phi(:,i) = phi(:,i)*modeFactor;
% end

fprintf('Calculating POD modes\n')
phi = u*V;
modeFactor = 1./ModeFactor;
phi = phi.*modeFactor';

phi = phi';
