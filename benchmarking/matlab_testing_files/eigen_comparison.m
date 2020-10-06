A = readmatrix('./matlab-mat.txt');
Sc = readmatrix('./matlab-membershipvec.txt'); %membership vector from the C implementation

tic
[V1, D1] = eig(A); 
toc
eigen_values = diag(D1);       % get all of the eigenvalues

[largest_eigen_value, index_of_largest] = max((eigen_values)); % find the largest eigenvalue
eig_v = V1(:,index_of_largest); % eigen vector that corresponds to the largest eigen value
eig_v * largest_eigen_value;
eig_v;

Sm = (eig_v > 0) - (eig_v < 0);

if (Sm(1) ~= Sc(1))
    Sm = Sm*Sc(1);

end
S_vecs_same = sum(Sc == Sm) == size(A,1)
differences = abs(sum(Sc == Sm) - size(A,1))
if (~S_vecs_same)
    sprintf('Not equal! Total differences:  %d', differences)
end
