function [ pars, pars_temp ] = pars_grey_ex1(t, xt_noise)
% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

    %% reponse vector and design matrix
    n = length(t); 
    
    h = [1; diff(t)];
    yt = cumsum(h .* xt_noise);                 % cumulative sum
    
    yt_trap = 1/2*yt(1:n-1) + 1/2*yt(2:n);      % discretization using trapezoid rule
    
    u_2 = (t(1:n-1).^2+t(2:n).^2)/2;            % independent varibale u_2(t)=t^2
    u_1 = (t(1:n-1)+t(2:n))/2;                  %                      u_1(t)=t
    init = ones(n-1,1);                         % constant value
    
    Theta = [yt_trap u_2 u_1 init];             % design matrix for regression
    
    X  = xt_noise(2:n);                         % response varibale 
    
    %% parameter estimation using regression formula 
    Xi = regress(X, Theta)';            % [a, b_2, b_1, c]
    
    %% initial value selection 
    eta = 1/(1-Xi(1))*(Xi(2)*t(1)^2+Xi(3)*t(1)+Xi(4)); 
    
    pars_temp = [Xi eta];               % [a, b_2, b_1, c, eta]
    
    pars = [Xi(1) 2*Xi(2) Xi(3) eta]; 
       
end