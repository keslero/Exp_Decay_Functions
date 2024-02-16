clc, clear, close all
% warning('off')

%% ADD PATHS
addpath('POD_SNP-25-50/')
addpath('POD_SNP-100-125/')
addpath('POD_SNP-125-150/')
addpath('POD_SNP-150-175/')
addpath('POD_SNP-175-200/')
addpath('POD_SNP-200-225/')
addpath('Functions/')
addpath('Figures/')

%% RUN 'trans_erg_t' INDUSTRIALLY
% -> choice (which energies to plot):
%        case 1 % 'r'     - Radial
%        case 2 % 't'     - Azimuthal
%        case 3 % 'z'     - Axial
%        case 4 % 'r-t'   - Transverse
%        case 5 % 'r-t-z' - Total
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% tm =  ...
%     [25 50
%     100 125
%     125 150
% 150 175
%     175 200
%     200 225];
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

choice = [1];

tm =  ...
    [100 125
    125 150
    150 175
    175 200];

start  = {[1 2 3] % Only  Mode   1,   2,   3
    [1 5 10]      % Only  Mode   1,   5,  10
    [1 1 1]       % First Modes  2,   3,   5
    [1 1 1]       % First Modes  5,  10, 100
    [1 2 3]};     % W/O   Modes  -,   1, 1-2
finish = {[1 2 3]
    [1 5 10]
    [2 3 5]
    [5 10 200]
    [200 200 200]};

cases = size(start,1);
for i = 1:cases
    np = cell2mat(start (i,:));
    nq = cell2mat(finish(i,:));

    fprintf('\n\nCase %d/%d\n\n',i,cases)
    main(tm,np,nq,choice)
end

fprintf('\nDone\n')