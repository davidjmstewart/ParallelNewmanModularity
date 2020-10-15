format compact
B_gc = readmatrix('./matlab-mat.txt'); % B_g matrix from the C implementation
Sc = readmatrix('./matlab-membershipvec.txt'); % B_g matrix from the C implementation
degrees = sum(A');

m = (sum(degrees))/2 ; % total number of edges in the network
n = size(A, 1);

K = degrees' * degrees;
B = (A - K/(2*m));
Q = 0;
ng = length(B);
B_gm = B - eye(ng, ng) .* sum(B'); % equation 6 in the paper


[V1, D1] = eig(B_gm); 

eigen_values = diag(D1);       % get all of the eigenvalues

[largest_eigen_value, index_of_largest] = max((eigen_values)); % find the largest eigenvalue
eig_v = V1(:,index_of_largest); % eigen vector that corresponds to the largest eigen value
Sm = (eig_v > 0) - (eig_v < 0);

deltaQ = Sm' * B_gm * Sm; % equation 5 in the paper

sum(sum(abs(B_gm - B_gc)))
sum(sum(abs(B_gm - B_gc))) < 0.001