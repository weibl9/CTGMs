function [ pars_hat ] = pars_integral_ex1(t, yt)
% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

% t: time point vector
% xt: time series data with noise 
% regression formua: 
    % x(t_k) = \beta(1) \int_{t_1}^{t_k} x(\tau)d\tau + 
    %          \beta(2) \int_{t_1}^{t_k} \tau d\tau + 
    %          \beta(3) \int_{t_1}^{t_k} 1 d\tau + 
    %          x(t_1)

h = diff(t);                                    % time difference
n = length(yt); 

cums_trap = cumsum(h.*(yt(1:n-1)+yt(2:n)))/2;   % trapezoid cusums 
lin_vec = cumsum(h.*(t(1:n-1)+t(2:n)))/2;       % time linear term
cons_vec = cumsum(h.*ones(n-1,1));              % constant term 
ini_vec = ones(n-1, 1);                         % initial value

Xt = [cums_trap lin_vec cons_vec ini_vec];      % design matrix for regression

pars_hat = regress(yt(2:n), Xt)';               % [beta, ini_val]
    
end