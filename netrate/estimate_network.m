function [A_hat, total_obj] = estimate_network(A, C, num_nodes, horizon, type_diffusion),

num_cascades = zeros(1,num_nodes);

% A_potential(j,i) = sum_{all cascades} ( t_i - t_j )
A_potential = sparse(zeros(size(A))); % A is nxn, n is the number of nodes

% A_bad(i,j) = sum_{all cascades} ( T - t_i )
A_bad = sparse(zeros(size(A)));
A_hat = sparse(zeros(size(A))); % Estimated pair-wise transmission rates
total_obj = 0;

for c=1:size(C, 1),
    idx = find(C(c,:)~=-1); % used nodes, idx has node indices
    [val, ord] = sort(C(c, idx)); % idx(ord(i)) will give node indices in increasing order of infection time
                                  % val(i) = C(c, idx(ord(i))) gives the infection time of node idx(ord(i))

    for i=2:length(val),
        num_cascades(idx(ord(i))) = num_cascades(idx(ord(i))) + 1;
        % Set A_potential(idx(ord(j)), idx(ord(i))) = sum of (t_i - t_j) for all t_j < t_i
        for j=1:i-1,
            if (strcmp(type_diffusion, 'exp'))
                A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+val(i)-val(j);
            elseif (strcmp(type_diffusion, 'pl') && (val(i)-val(j)) > 1)
                A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+log(val(i)-val(j));
            elseif (strcmp(type_diffusion, 'rayleigh'))
                A_potential(idx(ord(j)), idx(ord(i))) = A_potential(idx(ord(j)), idx(ord(i)))+0.5*(val(i)-val(j))^2;
            end
        end
    end
    
    for j=1:num_nodes,
        if isempty(find(idx==j)) % Node j is not infected in cascade c
            % Set A_bad(idx(ord(i)), j) = sum of (T - t_i) for all t_i in infected nodes in c
            for i=1:length(val), % val has the infection times of the infected nodes in cascade c
                if (strcmp(type_diffusion, 'exp'))
                    A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + (horizon-val(i));
                elseif (strcmp(type_diffusion, 'pl') && (horizon-val(i)) > 1)
                    A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + log(horizon-val(i));
                elseif (strcmp(type_diffusion, 'rayleigh'))
                    A_bad(idx(ord(i)), j) = A_bad(idx(ord(i)), j) + 0.5*(horizon-val(i))^2;
                end
            end
        end
    end
end

% we will have a convex program per column

for i=1:num_nodes,
    % Check how many cascades node i is a part of
    if (num_cascades(i)==0)
        A_hat(:,i) = 0; % Column i has transmission rates from all nodes to i
        continue;
    end
    
    % 
    cvx_begin
        variable a_hat(num_nodes);
        variable t_hat(num_cascades(i));
        obj = 0;
    
        a_hat(A_potential(:,i)==0) == 0; % If any node never infects i in any cascade, set its a_ij = 0
        
        for j=1:num_nodes,
            if (A_potential(j,i) > 0)
                obj = -a_hat(j)*(A_potential(j,i) + A_bad(j,i)) + obj; % TODO
            end
        end
        
        c_act = 1;
        for c=1:size(C, 1),
            idx = find(C(c,:)~=-1); % used nodes
            [val, ord] = sort(C(c, idx));
    
            idx_ord = idx(ord); % Node indices ordered by infection time
            cidx = find(idx_ord==i); % Node i is the cidx'th infected node
            
            if (~isempty(cidx) && cidx > 1) % Skip first node and non-infected nodes
                if (strcmp(type_diffusion, 'exp'))
                    t_hat(c_act) == sum(a_hat(idx_ord(1:cidx-1)));
                elseif (strcmp(type_diffusion, 'pl'))       
                    tdifs = 1./(val(cidx)-val(1:cidx-1));
                    indv = find(tdifs<1);
                    tdifs = tdifs(indv);
                    t_hat(c_act) <= (tdifs*a_hat(idx_ord(indv)));
                elseif (strcmp(type_diffusion, 'rayleigh'))
                    tdifs = (val(cidx)-val(1:cidx-1)); 
                    t_hat(c_act) <= (tdifs*a_hat(idx_ord(1:cidx-1)));
                end
                
                obj = obj + log(t_hat(c_act));
                
                c_act = c_act + 1;
            end
        end

        a_hat >= 0;
        
        maximize obj
    cvx_end
    
    total_obj = total_obj + obj;
    A_hat(:,i) = a_hat;
end
