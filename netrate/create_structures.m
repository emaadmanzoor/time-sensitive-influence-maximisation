function [num_cascades, A_potential, A_bad, A_hat] ...
        = create_structures(A, C, num_nodes, horizon, type_diffusion)
% TODO
% Create matrix S(c,i) = state of node i in cascade s
% Select states randomly from 1 to num_states

% TODO
% Create vector StartTime(c) = start time of cascade c

% TODO
% Create a num_cascades{s}, A_potential{s}, A_bad{s} and A_hat{s}
% for each state s. If a node is in A_potential{s}, it is in state s and has
% edges from all other nodes to it added to A_potential{s}(j, i)

%%
num_cascades = zeros(1,num_nodes);

% A_potential(j,i) = sum_{all cascades} ( ti2 - tj2 )
%A_potential = sparse(zeros(size(A))); % A is nxn, n is the number of nodes
A_potential = zeros(size(A)); % A is nxn, n is the number of nodes

% A_bad(i,j) = sum_{all cascades} ( T - ti2 )
%A_bad = sparse(zeros(size(A)));
%A_hat = sparse(zeros(size(A)));
A_bad = zeros(size(A));
A_hat = zeros(size(A));

disp 'Creating structures...'

for c=1:size(C, 1),
    idx = find(C(c,:)~=-1); % used nodes, idx has node indices
    [val, ord] = sort(C(c, idx)); % idx(ord(i)) will give node indices in increasing order of infection time
                                  % val(i) = C(c, idx(ord(i))) gives the infection time of node idx(ord(i))
    for i=2:length(val),
        i2 = idx(ord(i));
        num_cascades(i2) = num_cascades(i2) + 1;
        % Set A_potential(idx(ord(j)), idx(ord(i))) = sum of (ti2 - tj2) for all tj2 < ti2
        for j=1:i-1,
            j2 = idx(ord(j));
            %if (strcmp(type_diffusion, 'exp'))
                A_potential(j2, i2) = A_potential(j2, i2)+val(i)-val(j);
            %elseif (strcmp(type_diffusion, 'pl') && (val(i)-val(j)) > 1)
            %    A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+log(val(i)-val(j));
            %elseif (strcmp(type_diffusion, 'rayleigh'))
            %    A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+0.5*(val(i)-val(j))^2;
            %end
        end
    end
    
    for j=1:num_nodes,
        if isempty(find(idx==j, 1)) % Node j is not infected in cascade c
            % Set A_bad(idx(ord(i)), j) = sum of (T - ti2) for all ti2 in infected nodes in c
            for i=1:length(val), % val has the infection times of the infected nodes in cascade c
            %    if (strcmp(type_diffusion, 'exp'))
                    A_bad(i2, j) = A_bad(i2, j) + (horizon-val(i));
            %    elseif (strcmp(type_diffusion, 'pl') && (horizon-val(i)) > 1)
            %        A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + log(horizon-val(i));
            %    elseif (strcmp(type_diffusion, 'rayleigh'))
            %        A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + 0.5*(horizon-val(i))^2;
            %    end
            end
        end
    end
end

disp 'Done creating structures.'
