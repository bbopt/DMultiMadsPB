% As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
% Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary
% Computation 8(2): 173-195, 2000.
%
% Example T3.
%
% This file is part of a collection of problems developed for
% derivative-free multiobjective optimization in
% A. L. Cust?dio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
% Direct Multisearch for Multiobjective Optimization, 2010
% implemented in matlab.
%
% Written by the authors in June 1, 2019.
%
% x var : n = 30
% x <= 1.0
% x >= 0.0
%
function [F] = ZDT3(x);
    % params
    m = 30;

    % functions
    f(1) = x(1);

    gx = 1 + 9 / (m -1) * sum(x(2:m));
    h =  1- sqrt(f(1) / gx) - (f(1) / gx) * sin(10 *pi * f(1));

    f(2) = gx * h;

    F = f';
