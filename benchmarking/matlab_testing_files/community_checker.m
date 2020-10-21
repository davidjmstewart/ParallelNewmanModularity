close all;
clear all;
clc;
f = @(x)isequal(length(x), length(unique(x)));
A = readmatrix('./integer-matrix-matlab-mat.txt'); % B_g matrix from the C implementation

% A = recon_adjacency_matrix; % redefine adajency_matrix as A so that it matches Newman's notation
% A = jazzMat;
degrees = sum(A');

m = (sum(degrees))/2 ; % total number of edges in the network
n = size(A, 1);

K = degrees' * degrees;
B = (A - K/(2*m));
Q = 0;
tic
[communities, label, Q] = assignCommunity(B, 0, [1:length(A)]);
Q = Q/(4*m)
toc;
nodes(:,3) = num2cell(communities');
communities;



fid = fopen('matlab-communities.txt');
tline = fgetl(fid);
c = 0;
num_differences = 0;
while ischar(tline)
    c_communities = str2num(tline);
    matlab_communities = find(communities == c) - 1;
    
    if (length(setdiff(c_communities, matlab_communities)))
        disp('Not equal! Here is the difference')
        setdiff(c_communities, matlab_communities);
        num_differences = num_differences + (length(setdiff(c_communities, matlab_communities)));
    end
    tline = fgetl(fid);
    c = c + 1;
end


num_differences
fclose(fid);