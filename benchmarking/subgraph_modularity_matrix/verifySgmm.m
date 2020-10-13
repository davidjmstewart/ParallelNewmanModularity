ng = length(B);  % use this to construct the n_g x n_g subraph modualrity matrix (equation 6)
B_g = B - eye(ng, ng) .* sum(B'); % equation 6 in the paper

A = readmatrix('./matlab-mat.txt');
B_gc = readmatrix('./matlab-membershipvec.txt'); %membership vector from the C implementation
