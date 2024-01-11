%% Hausdorff Distance: Compute the Hausdorff distance between two point clouds.
%  Modified by Shenghua Ye based on:
%  Hassan R-Esfahlan (2010), Hausdorff Distance 
% (https://www.mathworks.com/matlabcentral/fileexchange/27905-hausdorff-distance), 
% MATLAB Central File Exchange. 

% Let A and B be subsets of a metric space (Z,dZ), 
% The Hausdorff distance between A and B, denoted by dH (A, B), is defined by:
% dH (A, B)=max{sup dz(a,B), sup dz(b,A)}, for all a in A, b in B,
% dH(A, B) = max(h(A, B),h(B, A)),  
% where h(A, B) = max(min(d(a, b))),  
% and d(a, b) is a L2 norm. 
% dist_H = hausdorff( A, B ) 
% A: First point sets. 
% B: Second point sets. 
% ** A and B may have different number of rows, but must have the same number of columns. ** 

%%
function [dH] = hausdorff( A, B, inv_cov) 
if(size(A,2) ~= size(B,2)) 
    fprintf( 'WARNING: dimensionality must be the same\n' ); 
    dist = []; 
    return; 
end
dH = max(compute_dist(A, B, inv_cov), compute_dist(B, A, inv_cov));

%% Compute distance
function[dist] = compute_dist(A, B, inv_cov) 
m = size(A, 1); 
n = size(B, 1); 
dim= size(A, 2); 
for k = 1:m 
    C = ones(n, 1) * A(k, :);
    D = diag( (C-B) * inv_cov * (C-B)' );
    D = sqrt(D);
%     D = (C-B) .* (C-B); 
%     D = sqrt(D * ones(dim,1)); 
    dist(k) = min(D); 
end
dist = max(dist);
