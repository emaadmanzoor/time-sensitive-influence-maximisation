function [A_hat, total_obj, pr, mae] = netrate_mixture(network, cascades, horizon, type_diffusion, num_nodes)
%%
min_tol = 1e-4;

pr = zeros(1,2);

disp 'Reading groundtruth...'
A_full = create_adj_matrix(network, num_nodes);

disp 'Reading cascades...'
C_full = create_cascades(cascades, num_nodes);

% Truncate graph and cascades
A = A_full(1:10, 1:10);
C = C_full(:, 1:10);

%%
disp 'Building data structures...'
[A_hat, total_obj] = estimate_network_mixture(A, C, num_nodes, horizon, type_diffusion);
%%

if exist(network),
    mae = mean(abs(A_hat(A~=0)-A(A~=0))./A(A~=0)); % mae
    pr(2) = sum(sum(A_hat>min_tol & A>min_tol))/sum(sum(A>min_tol)); % recall
    pr(1) = sum(sum(A>min_tol))/sum(sum(A_hat>min_tol)); % precision
else
    mae = [];
    pr = [];
    
end

save(['solution-', network], 'A_hat', 'mae', 'pr', total_obj);
