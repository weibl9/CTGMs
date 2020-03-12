function [ Aim, Condim, Xt_im_fits ] = odes_origin(t, Xt, tspan, solver)
% Copyright 2020, All Rights Reserved
% Code by Baolei Wei
% For paper, "On novel framework for continuous-time grey models: 
%                an integral matching perspective"
% by Baolei Wei, Naiming Xie

%Integral matching for original series
%   xt = [xt1(t_1:t_n), xt2(t_1:t_n)] 
%   tspan: time index including [ttrain; ttest]
%   model formula: d/dt xt = A*xt

[nobs, nsta] = size(Xt);                    % size of Xt

%% discrete integral transformation
h = diff(t); 
Yt = [h h] .* cumsum(Xt(1:nobs-1,:)+Xt(2:nobs,:))/2;

%% parameter estimation: two column with two regression 
Const = ones(nobs-1,1);

Pi = zeros(nsta, nsta+1);
for ind_stat = 1:nsta
    pars = regress(Xt(2:nobs,ind_stat), [Yt, Const]);
    Pi(ind_stat,:) = pars';
end

Aim = Pi(:,1:nsta); Condim = Pi(:,1+nsta);

if solver == 1
    %% method II: the expression of the close-form solution (FAST)
    % obtain the closed-form solution 
    % syms x1(t) x2(t)
    % syms a11 a12 a21 a22 cond1 cond2
    % x = [x1; x2];
    % Aim = [a11 a12; a21 a22]; Condim = [cond1; cond2];
    % odes = diff(x) == Aim*x; cond = x(0) == Condim;
    % [x1Sol(t), x2Sol(t)] = dsolve(odes, cond);
    % x1Sol(t)
    % x2Sol(t)
    
    a11 = Aim(1,1); 
    a12 = Aim(1,2); 
    a21 = Aim(2,1); 
    a22 = Aim(2,2); 
    cond1 = Condim(1); 
    cond2 = Condim(2);

    nobs = length(tspan);
    Xt_im_fits = zeros(nobs, nsta);
    for iter_t = 1:nobs
        t = tspan(iter_t);
        Xt_im_fits(iter_t,:) = ...
        [- (exp((t*(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*(a22/a21 - (a11/2 + a22/2 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)/2)/a21)*(2*a21*cond1 - a11*cond2 + a22*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) - (exp((t*(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*(a22/a21 - (a11/2 + a22/2 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)/2)/a21)*(a11*cond2 - 2*a21*cond1 - a22*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)), ...
         (exp((t*(a11 + a22 - (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*(a11*cond2 - 2*a21*cond1 - a22*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) + (exp((t*(a11 + a22 + (a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/2)*(2*a21*cond1 - a11*cond2 + a22*cond2 + cond2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)))/(2*(a11^2 - 2*a11*a22 + a22^2 + 4*a12*a21)^(1/2)) ...
        ];
    end
else
    %% method I: close-form soluion to ODEs d/dt x = A*x
    syms x1(t) x2(t)
    x = [x1; x2]; 
    odes = diff(x) == Aim*x; cond = x(0) == Condim; 
    [x1Sol(t), x2Sol(t)] = dsolve(odes, cond);

    %% fitting and forecasting values of orgin series
    Xt_im_fits = [eval(subs(x1Sol(t), t, tspan)), ...
                eval(subs(x2Sol(t), t, tspan))];
end

end

