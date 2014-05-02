%%
network = 'kronecker-core-periphery-n1024-h10-r0_01-0_25-network.txt';
cascades = 'kronecker-core-periphery-n1024-h10-r0_01-0_25-1000-cascades.txt';

horizon = 10;
type_diffusion = 'exp';
num_nodes = 10;
%%

[A_hat, total_obj, pr, mae] = netrate_mixture(network, cascades, horizon, type_diffusion, num_nodes)