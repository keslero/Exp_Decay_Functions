clc, clear, close all

% Choose format to change to
fig_format = '.png';

start  = [1 1 1
    1 1 1
    1 2 5
    2 2 2
    2 2 2
    1 2 5];
finish = [1 2 5
    10 100 200
    1 2 5
    2 5 10
    2 10 200
    200 200 200];

cases = size(start,1);

for i = 1:cases
    np = start(i,:);
    nq = finish(i,:);

    mode_combos = length(np);
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
    fig_fig = strcat(figName,'.fig');
    fig_ID = open(fig_fig);
    fig_ID.WindowState = 'maximized';
    fig_png = strcat('Figures/',figName,fig_format);
    saveas(fig_ID,fig_png)
    close(fig_ID)
end