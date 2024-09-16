#include "Graph.h"

#include "assert.h"
#include <string>
#include <sstream>
#include <vector>

void ReadClusterOutput(const std::string &, std::vector<std::size_t> &, std::size_t &);
void ReadNodeNumber(const std::string &, std::vector<std::size_t> &);
void ConvertLocalClusterNum(const std::vector<std::size_t> &, const std::size_t, const std::vector<std::size_t> &, std::vector<std::size_t> &, std::size_t &);
void RecordResults(std::string &, std::vector<std::size_t>&, const int, const int);

int main(int argc, char **argv) {
  // Check the input format
	if (argc != 6) {
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage %s orig.gml output_dir output_tag num_comps obj_name\n", argv[0]);
		exit(1);
	}

  // Save inputs to local variables
  std::string graph_file(argv[1]);
	std::string out_dir(argv[2]);
	std::string out_tag(argv[3]);
  int num_comps = atoi(argv[4]);
  std::string obj_name(argv[5]);

  Graph my_graph;
  my_graph.readFile(graph_file);
	
  std::vector<std::size_t> nn_comp, local_clust_num, cluster_num(my_graph.getNumNodes());
  std::size_t cluster_count = 0, local_cluster_count, num_comps_split = 0;

  std::string comp_filename = out_dir + out_tag + "_comp0"; 
  std::string out_filename = comp_filename + ".out";
  std::string nn_filename = comp_filename + ".nn";

  ReadNodeNumber(nn_filename, nn_comp); // Read .nn file for singletons
  if(nn_comp.size() > 0) {
    local_clust_num.resize(nn_comp.size());
    for(std::size_t i = 0; i < nn_comp.size(); ++i)
      local_clust_num[i] = i + 1;

    ConvertLocalClusterNum(local_clust_num, nn_comp.size(), nn_comp, cluster_num, cluster_count);
  }

  // Loop through all components
  for(std::size_t comp = 1; comp <= num_comps; ++comp) {
    comp_filename = out_dir + out_tag + "_comp" + std::to_string(comp); 
    out_filename = comp_filename + "_" + obj_name + ".out";
    nn_filename = comp_filename + ".nn";

    // Read .out file
    ReadClusterOutput(out_filename, local_clust_num, local_cluster_count);

    if(local_cluster_count > 1) {
      printf("Comp %lu split\n", comp);
      ++num_comps_split;
    }
    // Read .nn file
    ReadNodeNumber(nn_filename, nn_comp);

    // Set .out results based on .nn data
    ConvertLocalClusterNum(local_clust_num, local_cluster_count, nn_comp, cluster_num, cluster_count);        
  }
  printf("Number of components split: %lu\n\n\n", num_comps_split);

  // Record Results
  std::string output_file = out_dir + out_tag + "_" + obj_name +".out";
  RecordResults(output_file, cluster_num, cluster_count, my_graph.getNumEdges());

  return 0;
}

void ReadClusterOutput(const std::string &clust_file,
                       std::vector<std::size_t> &clust_num,
                       std::size_t &cluster_count) {
  FILE* input;
	char tmp_str[50];
  std::size_t i = 0, num_nodes, num_edges;
  int tmp_num;
	tmp_str[0] = '\0';
	
  if((input = fopen(clust_file.c_str(), "r")) == NULL) {
		printf("ERROR - Cound not open clust.out file (%s)\n.", clust_file.c_str());
		exit(1);
	}

  fscanf(input, "%lu", &num_nodes); // number of nodes
  fscanf(input, "%s", tmp_str); // string
  fscanf(input, "%lu", &cluster_count); // number of clusters    
  fscanf(input, "%s", tmp_str); // string
  fscanf(input, "%lu", &num_edges); // number of edges
  fscanf(input, "%s", tmp_str); // strings
    
  clust_num.resize(num_nodes);

  while(1) {
    // Check if end-of-file was reached
		if (feof(input)) {
			break; // Break while loop
		}

		if(fscanf(input, "%i", &tmp_num) > 0) {
      if(tmp_num == -1)
        clust_num[i++] = ++cluster_count;
      else
        clust_num[i++] = (std::size_t)tmp_num;
    }
	}

	fclose(input);    
}

void ReadNodeNumber(const std::string &nn_file,
                    std::vector<std::size_t> &node_num) {
  FILE* input;
	std::size_t node_id, num_nodes = 0;    
		
	if((input = fopen(nn_file.c_str(), "r")) == NULL) {
		printf("ERROR - Cound not open nn file.");
		exit(1);
	}

  while(1) {
    // Check if end-of-file was reached
    if (feof(input))
      break; // Break while loop
            
    if(fscanf(input, "%lu", &node_id) > 0)
      ++num_nodes;
	}

  rewind(input);
  node_num.resize(num_nodes);

  std::size_t i = 0;
  while(1) {
    // Check if end-of-file was reached
    if (feof(input))
      break; // Break while loop
            
    if(fscanf(input, "%lu", &node_id) > 0)
      node_num[i++] = node_id;
	}

	fclose(input);    
}

void ConvertLocalClusterNum(const std::vector<std::size_t> &local_cluster_num,
                            const std::size_t local_cluster_count,
                            const std::vector<std::size_t> &nn,
                            std::vector<std::size_t> &cluster_num,
                            std::size_t &cluster_count) {
  assert(local_cluster_num.size() == nn.size());

  // Loop through all nodes in local_cluster_num
  for(std::size_t i = 0; i < local_cluster_num.size(); ++i)
    cluster_num[nn[i]] = local_cluster_num[i] + cluster_count;
    
  cluster_count += local_cluster_count;
}

void RecordResults(std::string &output_file,
                   std::vector<std::size_t>& clust,
                   const int num_clusts,
                   const int num_edges) {
  FILE* output;
  if((output = fopen(output_file.c_str(), "w")) == NULL) {
    printf("ERROR - Could not open file (%s)\n", output_file.c_str());        
  }

  fprintf(output, "%lu nodes %i clusters %i edges\n", clust.size(), num_clusts, num_edges);
  for(std::size_t i = 0; i < clust.size(); ++i)
    fprintf(output, "%d ", clust[i]);
  fclose(output);
}
