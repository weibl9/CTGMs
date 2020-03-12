% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

clc; clear; close
addpath('./utils');
addpath('./results');

%% set initial parameters
nrep = 1000;                        % number of replications
hSet = [0.25 0.10 0.05];            % time interval (sample size)
snrSet = [2.5 3.5 5.0];             % signal-to-noise (snr) ratio

ini = 1.2;                          % initial value of ODE
str = [0.15, 0.20, -0.25];          % parameter vector of ODE

mean4pars = zeros(9, 10);           % sample means of estimated parameters
std4pars = zeros(9, 10);            % standrad deviation of estimated parameters 

%% main loop 
for iter_h=1:length(hSet)
    h = hSet(iter_h);               % each sample size     
    for iter_snr = 1:length(snrSet)        
        snr = snrSet(iter_snr);     % each snr ratio         
        %% nosie-free time series
        tim_train = 0.0:h:5.0;               % time instants for training
        nobs_train = length(tim_train); 
        nobs_test = 10;                      % number of test samples 
        tim_test = 5.0 + (1:nobs_test)*h;    % time instants for testing 
        
        tspan = [tim_train tim_test]; 
        [~, xt] = ode45(@(t,x)ode_im(t,x,str),tspan,ini); % numerical solution

        %% 1000(nrep) replications in each sample size and snr combination
        pars_mat_im = zeros(nrep,4); pars_mat_gm = zeros(nrep,4);   % save estimated parameters
        mape_mat_im = zeros(nrep,4); mape_mat_gm = zeros(nrep,4);   % save mapes in each replication
        for iter_rep=1:nrep
            %% add noise to noise-free time series
            rng(iter_rep);                              % repeatable random numbers 
            sigma = 1/(snr^2)*std(xt,1);                % snr ratio: SNR=sqrt(1/(1/4))=2
            xt_noise = xt(1:nobs_train)+normrnd(0,sigma,[nobs_train, 1]);

            %% integral matching 
            pars_im = pars_integral_ex1(tim_train',xt_noise);        % parameter estimation
            [~, xt_im] = ode45(@(t,x)ode_im(t,x,pars_im(1:3)),tspan,pars_im(4));  % fits and forecasts

            pars_mat_im(iter_rep,:) = pars_im;         % estimated parameters in each replication 
            ape = abs( (xt-xt_im)./xt )*100;           % fitting and forecasting in each replication
            mape_mat_im(iter_rep,:) = [ 
                mean(ape(1:nobs_train)), ...           % fitting
                mean(ape(nobs_train+(1:2))), ...       % 2-step ahead 
                mean(ape(nobs_train+(1:5))), ...       % 5-step ahead 
                mean(ape(nobs_train+(1:nobs_test))) ]; % nobs_test-step ahead     

            %% grey modelling 
            [pars_gm,pars_temp] = pars_grey_ex1(tim_train',xt_noise); % parameter estimation
            [~, yt_gm] = ode45(@(t,x)ode_gm(t,x,pars_temp(1:4)),tspan,pars_temp(5));
            hspan = (diff(tspan))';
            xt_gm = [ pars_temp(5); 
                     (yt_gm(2:length(yt_gm)) - yt_gm(1:length(yt_gm)-1))./hspan]; % fits and forecasts

            pars_mat_gm(iter_rep,:) = pars_gm;         % estimated parameters in each replication 
            ape = abs( (xt-xt_gm)./xt )*100;           % fitting and forecasting in each replication
            mape_mat_gm(iter_rep,:) = [ 
                mean(ape(1:nobs_train)), ...           % fitting
                mean(ape(nobs_train+(1:2))), ...       % 2-step ahead
                mean(ape(nobs_train+(1:5))), ...       % 5-step ahead 
                mean(ape(nobs_train+(1:nobs_test))) ]; % nobs_test-step ahead

        end
        %% sample mean and standard derivation of estimated parameters
        mean4pars(length(hSet)*(iter_h-1)+iter_snr,:) = [5/h+1, snr, mean(pars_mat_gm), mean(pars_mat_im)];
        std4pars(length(hSet)*(iter_h-1)+iter_snr,:) = [5/h+1, snr, std(pars_mat_gm), std(pars_mat_im)];
        
        %% mape values for FIGURE 4!
        name = sprintf('mape_%d_%d_gm.csv',10*snr,5/h+1); 
        csvwrite(['results/',name], mape_mat_gm)
        name = sprintf('mape_%d_%d_im.csv',10*snr,5/h+1); 
        csvwrite(['results/',name], mape_mat_im)        
    end
end

%% Table 2!
format short g
round(mean4pars,3)
round(std4pars, 3)












