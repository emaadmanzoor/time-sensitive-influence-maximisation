function [A_hat, total_obj] = estimate_network(A, C, num_nodes, horizon, type_diffusion, ...
                                               num_cascades, A_potential, A_bad, A_hat)
%%
total_obj = 0;

% E-step

% Create matrix SCurrent, SPrev to track the number of changed states

% Evaluate tau(s), sigma(s) for the priors from StartTimes
%
% for each state s:
%   find the tau(s) = MLE estimate of tau from StartTimes{s}
%   find the sigma(s) = MLE estimate of sigma from StartTimes{s}

% Evaluate f(\alpha_ij | state = s) for each s \in states
% Results will be stored in A_hat itself
%
% for each state s
%   Use A_potential{s}, A_bad{s}
%   for i in nodes
%     -- convex opt --
%   end
%   Write to A_hat{s}
% end

% M-step

% Reassign states for all nodes
% 
% for c in cascades:
%   for i in nodes of cascade c:
%     likelihoods(s) = compute probability of s_ic given C(c,j)'s and alpha_hat
%     priors(s) = compute prior of s_ic given tau(s), sigma(s)
%     posteriors(s) = likelihoods * priors
%     S(c,i) = argmax(posteriors)
%   end
% end

% we will have a convex program per column

for i=1:num_nodes,
    % Check how many cascades node i is a part of
    if (num_cascades(i)==0)
        A_hat(:,i) = 0; % Column i has transmission rates from all nodes to i
        continue;
    end
    
    fprintf('Solving convex program for node %d\n', i);
    
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
