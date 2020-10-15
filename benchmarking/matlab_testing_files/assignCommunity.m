function [labelled_vertices, label, deltaQ ] = assignCommunity(B, label, globalVertices)
%assignCommunity Assigns all vertices in B to a community by using Newman's
%modularity algorithm
    THRESHOLD = 0.0001;
    OPTIMISATION_THRESHOLD = 1;
    
    ng = length(B);  % use this to construct the n_g x n_g subraph modualrity matrix (equation 6)
    B_g = B - eye(ng, ng) .* sum(B'); % equation 6 in the paper
    [V1, D1] = eig(B_g);           % calculate the eigenvalues and eigenvectors of this matrix
    eigen_values = diag(D1);       % get all of the eigenvalues
    [largest_eigen_value, index_of_largest] = max(eigen_values); % find the largest eigenvalue
    eig_v = V1(:,index_of_largest); % eigen vector that corresponds to the largest eigen value
    s = (eig_v >= 0) - (eig_v < 0); % membership vector
   
    deltaQ = s' * B_g * s; % equation 5 in the paper
    deltaQ_updated = false;
    
    
    % small deltaQ indicates it is not a good partition for the nodes
    if (deltaQ <= THRESHOLD) % check to see if it's "close enough" to 0
              
        if (deltaQ_updated == false)
            deltaQ = 0; % we are not making the split so don't change the modularity
        end
        
        labelled_vertices = zeros(1,ng); % preallocate our array we will fill with our label
        labelled_vertices(:) = label;
        label = label + 1;
    else
        left_indices  = find(s == -1);
        right_indices = find(s == 1);

        left_B  = B(left_indices,left_indices);
        right_B = B(right_indices,right_indices);

        %recurse for the first group, try and split them up using the same appraoch
        [left, label, deltaQ1] =  assignCommunity(left_B, label, globalVertices(left_indices));
        
        %recurse for the second group, try and split them up using the same appraoch
        [right, label, deltaQ2] = assignCommunity(right_B, label, globalVertices(right_indices));
        
        deltaQ = deltaQ + deltaQ1 + deltaQ2;
        labelled_vertices = zeros(1, length(left) + length(right));
        labelled_vertices(left_indices) = left;
        labelled_vertices(right_indices) = right;
    end
end