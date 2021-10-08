%% Figures for Two Cases
clc; 
clear;
close all;
addpath('./utils')
addpath('./results')

%% load data for multi boxplot
% Fugure (a) for the first component 
for f = 1:3 
    switch f
         case 1
            % case 1
            mape3gm = csvread('x1_mape_3_21_gm.csv'); mape3im = csvread('x1_mape_3_21_im.csv');
            mape4gm = csvread('x1_mape_4_21_gm.csv'); mape4im = csvread('x1_mape_4_21_im.csv');
            mape5gm = csvread('x1_mape_5_21_gm.csv'); mape5im = csvread('x1_mape_5_21_im.csv');
         case 2
            % case 2
            mape3gm = csvread('x1_mape_3_51_gm.csv'); mape3im = csvread('x1_mape_3_51_im.csv');
            mape4gm = csvread('x1_mape_4_51_gm.csv'); mape4im = csvread('x1_mape_4_51_im.csv');
            mape5gm = csvread('x1_mape_5_51_gm.csv'); mape5im = csvread('x1_mape_5_51_im.csv');
         case 3
            % case 3
            mape3gm = csvread('x1_mape_3_101_gm.csv'); mape3im = csvread('x1_mape_3_101_im.csv');
            mape4gm = csvread('x1_mape_4_101_gm.csv'); mape4im = csvread('x1_mape_4_101_im.csv');
            mape5gm = csvread('x1_mape_5_101_gm.csv'); mape5im = csvread('x1_mape_5_101_im.csv');
    end

    for loop = 1:4
        % Create example data 
        A = log10([mape3gm(:,loop) mape4gm(:,loop) mape5gm(:,loop)]); 
        C = log10([mape3im(:,loop) mape4im(:,loop) mape5im(:,loop)]); 

        % prepare data 
        data = cell(3,3); 
        for ii = 1:size(data,1) 
            Ac{ii} = A(:,ii); 
            Cc{ii} = C(:,ii); 
        end 
        data = vertcat(Ac, Cc);

        xlab = {'2.5','3.5','5.0'}; 

        col = [0, 255, 0, 175; 
              255, 255,0, 175]/255;

        %%    
        subplot(3,4,4*(f-1)+loop)
        multiple_boxplot(data', xlab, {'grey modelling','integral matching'}, col')

        grid on
        grid minor
        ylim([-3.5 3.5])
        yticks([-3:1:3]); % ytickformat('%.2f')
        set(gca, 'fontsize', 13) % 'FontName', 'times', 
    end

end 
set(gcf,'Position',[50 90 1400 760])
pause; 
h = figure(1);
set(h,'PaperSize',[15 10]); 
print(h,'fig3a','-dpdf','-painters') % set the paper size then print/save it 

% Figure (b) for the second component
clf
for f = 1:3 
    switch f
         case 1
            % case 1
            mape3gm = csvread('x2_mape_3_21_gm.csv'); mape3im = csvread('x2_mape_3_21_im.csv');
            mape4gm = csvread('x2_mape_4_21_gm.csv'); mape4im = csvread('x2_mape_4_21_im.csv');
            mape5gm = csvread('x2_mape_5_21_gm.csv'); mape5im = csvread('x2_mape_5_21_im.csv');
         case 2
            % case 2
            mape3gm = csvread('x2_mape_3_51_gm.csv'); mape3im = csvread('x2_mape_3_51_im.csv');
            mape4gm = csvread('x2_mape_4_51_gm.csv'); mape4im = csvread('x2_mape_4_51_im.csv');
            mape5gm = csvread('x2_mape_5_51_gm.csv'); mape5im = csvread('x2_mape_5_51_im.csv');
         case 3
            % case 3
            mape3gm = csvread('x2_mape_3_101_gm.csv'); mape3im = csvread('x2_mape_3_101_im.csv');
            mape4gm = csvread('x2_mape_4_101_gm.csv'); mape4im = csvread('x2_mape_4_101_im.csv');
            mape5gm = csvread('x2_mape_5_101_gm.csv'); mape5im = csvread('x2_mape_5_101_im.csv');
    end

    for loop = 1:4
        % Create example data 
        A = log10([mape3gm(:,loop) mape4gm(:,loop) mape5gm(:,loop)]); 
        C = log10([mape3im(:,loop) mape4im(:,loop) mape5im(:,loop)]); 

        % prepare data 
        data = cell(3,3); 
        for ii = 1:size(data,1) 
            Ac{ii} = A(:,ii); 
            Cc{ii} = C(:,ii); 
        end 
        data = vertcat(Ac, Cc);

        xlab = {'2.5','3.5','5.0'}; 

        col = [0, 255, 0, 175; 
              255, 255,0, 175]/255;

        %%    
        subplot(3,4,4*(f-1)+loop)
        multiple_boxplot(data', xlab, {'grey modelling','integral matching'}, col')

        grid on
        grid minor
        ylim([-3.5 3.5])
        yticks([-3:1:3]); % ytickformat('%.2f')
        set(gca, 'fontsize', 13) % 'FontName', 'times', 
    end

end 
set(gcf,'Position',[50 90 1400 760])
pause; 
h = figure(1);
set(h,'PaperSize',[15 10]); 
print(h,'fig3b','-dpdf','-painters') % set the paper size then print/save it 



 






