close all;
clear all;
clc;

A = readmatrix('./integer-matrix-matlab-mat.txt'); % B_g matrix from the C implementation

% A = recon_adjacency_matrix; % redefine adajency_matrix as A so that it matches Newman's notation
% A = jazzMat;
degrees = sum(A');

m = (sum(degrees))/2 ; % total number of edges in the network
n = size(A, 1);

K = degrees' * degrees;
B = (A - K/(2*m));
Q = 0;

[communities, label, Q] = assignCommunity(B, 0, [1:length(A)]);
Q = Q/(4*m);
nodes(:,3) = num2cell(communities');
communities