% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

clc; clear; close
addpath('./utils');
addpath('./results');

%% closed-form solutions for ODE  
syms x1(t) x2(t)    % time variable 
x = [x1; x2];       % state variables 

% true parameters 
A = [-0.25 0.70; 0.75 -0.25];
cond = x(0) == [1.20; 0.35];

% system of odes with closed-form solution for Original ODE
odes = diff(x) == A*x; 
[x1Sol(t), x2Sol(t)] = dsolve(odes, cond);

%% simulation experiments 

h_set = [0.25 0.10 0.05]; %  
snr_set = [2.5 3.5 5.0]; %  
pars_sample_means = zeros(9, 14); % length(snr_set)*length(snr_set) rows
pars_sample_stds = zeros(9, 14);  % length(snr_set)*length(snr_set) rows

for iter_h = 1:length(h_set)    
    %% generating time series without noise
    h = h_set(iter_h);
    
    ttrain = (0:h:5.0)';
    ttest = 5.0 + (1:10)'*h;
    tspan = [ttrain; ttest];

    Xt_true = [eval(subs(x1Sol(t), t, tspan)),...
                eval(subs(x2Sol(t), t, tspan))];
    Xt_train = Xt_true(1:length(ttrain),:);
   
    %% adding noise to the true time series and the results 
    for iter_snr = 1:length(snr_set)        
        snr = snr_set(iter_snr);    
        nrep = 1000;                              % the number of replications 

        pars_gm = zeros(nrep,6); pars_im = zeros(nrep,6); 
        x1_mape_gm = zeros(nrep,4); x1_mape_im = zeros(nrep,4); 
        x2_mape_gm = zeros(nrep,4); x2_mape_im = zeros(nrep,4);    

        for iter_r = 1:nrep

            [ntrain, nsta] = size(Xt_train);
            Xt_noise = zeros(ntrain, nsta);
            Xt_train_std = std(Xt_train, 1);    

            % adding noise to the true time series: each component 
            for i = 1:nsta  
                rng(iter_r);
                sigma = 1/(snr^2)*Xt_train_std(i); % snr = sqrt(1/(1/4))=2
                Xt_noise(:,i) = Xt_train(:,i) + normrnd(0, sigma, [ntrain, 1]);
            end

            % grey modelling for vector series with results 
            [Agm, Cgm, Condgm, Xt_gm_fits] = odes_cusum(ttrain, Xt_noise, tspan, 1);
            pars_gm(iter_r,:) = [reshape(Agm',1,[]), Condgm'];
            
            % integral matching for vector series with results
            [Aim, Condim, Xt_im_fits] = odes_origin(ttrain, Xt_noise, tspan, 1);
            pars_im(iter_r,:) = [reshape(Aim',1,[]), Condim'];
            
            %  summary of results obtained from gm(grey modelling) and im(integral macting) 
            err_gm = abs(Xt_gm_fits - Xt_true)./Xt_true*100;
            x1_mape_gm(iter_r,:) = ...
                [mean(err_gm(1:ntrain,1)), mean(err_gm(ntrain+(1:2),1)), ...
                 mean(err_gm(ntrain+(1:5),1)), mean(err_gm(ntrain+(1:10),1))];
            x2_mape_gm(iter_r,:) = ...
                [mean(err_gm(1:ntrain,2)), mean(err_gm(ntrain+(1:2),2)), ...
                 mean(err_gm(ntrain+(1:5),2)), mean(err_gm(ntrain+(1:10),2))];

            err_im = abs(Xt_im_fits - Xt_true)./Xt_true*100;
            x1_mape_im(iter_r,:) = [mean(err_im(1:ntrain,1)), ... % fitting error
                mean(err_im(ntrain+(1:2),1)), ...                 % 2-step ahead 
                 mean(err_im(ntrain+(1:5),1)), ...                % 5-step ahead 
                 mean(err_im(ntrain+(1:10),1))];                  % 10-step ahead 
            x2_mape_im(iter_r,:) = [mean(err_im(1:ntrain,2)), ...
                mean(err_im(ntrain+(1:2),2)), ...
                 mean(err_im(ntrain+(1:5),2)), ...
                 mean(err_im(ntrain+(1:10),2))];
        end
        
        %% sample means and standard derivation of estimated parameters
        pars_sample_means(length(h_set)*(iter_h-1)+iter_snr,:) = [5/h+1, snr, mean(pars_gm), mean(pars_im)];
        pars_sample_stds(length(h_set)*(iter_h-1)+iter_snr,:) = [5/h+1, snr, std(pars_gm), std(pars_im)];
        
        %% mape values for FIGURE 5 and 6! 
        gm_name = sprintf('x1_mape_%d_%d_gm.csv', 10*snr, 5/h+1); 
        csvwrite(['results/',gm_name], x1_mape_gm)
        gm_name = sprintf('x2_mape_%d_%d_gm.csv', 10*snr, 5/h+1); 
        csvwrite(['results/',gm_name], x2_mape_gm)
        
        im_name = sprintf('x1_mape_%d_%d_im.csv', 10*snr, 5/h+1); 
        csvwrite(['results/',im_name], x1_mape_im)
        im_name = sprintf('x2_mape_%d_%d_im.csv', 10*snr, 5/h+1); 
        csvwrite(['results/',im_name], x2_mape_im)
    end    
end

%% Table 4!
format short g
round(pars_sample_means,3)
round(pars_sample_stds,3)




