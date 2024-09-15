using NetworkInference
using LightGraphs

algorithm = PIDCNetworkInference()
dataset_name = string(ARGS[1])

genes = get_nodes(dataset_name);

network = InferredNetwork(algorithm, genes);

write_network_file(string(ARGS[2]), network);
