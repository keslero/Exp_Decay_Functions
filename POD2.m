function [phi, at, lam, Xmean, nbasis, tenergy]=POD2(X,cenergy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the POD basis vectors 
%
% Function inputs:
% X       :  (n x nsnap) snapshot matrix
%            n = number of states in full space
%            nsnap = number of snapshots
% cenergy :  how much energy of the ensemble you want 
%            to capture, i.e. cenergy = 99.9 (percent)
% 
% Function outputs:
% phi     :   (n x nsnap) matrix containing POD basis vectors
% lam     :   (nsnap x 1) vector containing POD eigenvalues
% Xmean   :   (n x 1) the mean of the ensemble X
% nbasis  :   number of POD basis vectors you should use to capture cenergy
%             [phi_u, at_u, lam_u, u_mean, nbasis_u, tenergy_u]
%             =POD2(u(:,1:500),95); (percent) of the ensemble
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the mean of the snapshots
nsnap = size(X,2);
Xmean = sum(X,2)/nsnap;
% Obtain new snapshot ensemble with zero mean
for i=1:nsnap
    X1(:,i) = X(:,i)-Xmean;
end
%   METHOD OF SNAPSHOTS
% calculate the empirical correlation matrix C
%C = X1'*X1/nsnap;
% without substruction the mean field:
C = X'*X/nsnap;
% Calculate the POD basis
[at,evalueC] = eig(C);
%phi = X1 * evectorC;
%without substructing the mean field:
phi = X * at;
% Normalize (or not) the POD basis
for i=1:nsnap
    %phi(:,i) = phi(:,i)/norm(phi(:,i),2);
    phi(:,i) = phi(:,i);
end
% return the POD eigenvalues in a vector
lam = diag(evalueC);

% Rearrange POD eigenvalues, vectors in descending order.
% Note that the correlation matrix C is symmetric, so SVD and EIG
% will give the same evectorC and evalueC but they are already in
% descending order and hence we don't need to rearrange evectorC and evalueC
% if SVD is used
lam = rot90(lam,2);

c= real(lam);
save c
%phi = fliplr(phi);
%%%%  Find the number of POD basis capturing cenergy (percent) of energy
% compute the total energy
tenergy = sum(lam);
energy = 0.;
nbasis = 0;
i = 1;
while (((energy/tenergy)*100) < cenergy)
    energy = energy + lam(i);
    i = i+1;
end
nbasis = i;
        


