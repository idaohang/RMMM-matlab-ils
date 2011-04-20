function [U,W,R] = transform (A)
% This function calculates the U_j, W_j and R_j in the CC algorithm
% using QR factorization
%
global CC_m;
[Q,R] = qr(A);
Q_t = Q';
R = R(1:3,:);
U = Q_t(1:3,:);
W = Q_t(4:CC_m-1,:);
