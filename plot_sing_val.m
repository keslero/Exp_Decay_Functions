function plot_sing_val(S)
%% Plot Singular Values
if size(S,2) >1
    sigma = diag(S);
else
    sigma = S;
end
sigma = real(sigma); % Singular values (vector)
epsilon = 1e-2; % Lower limit of energy to plot [%]
N = length(sigma);
energy = 100*diag(sigma)/sum(diag(sigma)); % Relative energy of each
% singular value (out of 100%)

yyaxis left
semilogy(energy,'ko','Linewidth',2), hold on
grid on, ylim([epsilon 100])
xlabel('\sigma_n'), ylabel('E [%]')

tot_eng = zeros(N,1);
for i = 1:N
    tot_eng(i) = sum(energy(1:i));
end
erg_levels = [95;97.5;99];
n = length(erg_levels);
erg_sig = zeros(n,1);
for i = 1:n
    erg_sig(i) = find(tot_eng >= erg_levels(i),1);
end

yyaxis right
plot(1:N,tot_eng,'k-','LineWidth',1), hold on
plot(1:N,erg_levels.*ones(1,N),'r--'), hold on
for i = 1:n
    plot(erg_sig(i)*[1 1],[0 erg_levels(i)],'b--'), hold on
end
x_max = max(find(energy>epsilon))+2;
grid on, xlim([1 x_max]), ylim([0 100])
xlabel('\sigma_n'), ylabel('Total Energy [%]')
title({'Singular Value Energy Content';'&';'Accumulated Energy Content [%]'})

ax = gca;
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0 0 0];
end