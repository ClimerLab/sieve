#include "Graph.h"
#include <cstdio>
#include <cstdlib>
#include <string>

int main(int argc, char **argv) {	
	// Check inputs
	if (argc < 4) {
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage %s <input.gml> <output_dir> <output_tag>\n", argv[0]);
		exit(1);
	}

	// Read in user input to local variables
	std::string input_file(argv[1]);
	std::string output_dir(argv[2]);
	std::string output_tag(argv[3]);

	Graph my_graph;
		
	// Read in graph
	my_graph.readFile(input_file);	

	FILE* output;
	std::string num_comp_file = output_dir + "NUM_COMP.txt";
	if((output = fopen(num_comp_file.c_str(), "w")) == NULL) {
		printf("ERROR in SeperateGraphIntoComponents - Could not read in graph.\n");
		exit(1);
	}
	
	// Run BFS on graph and record a .gml and .nn file for each component
	fprintf(output, "NUM_COMP=%lu\n", my_graph.seperateGraphIntoComponents(output_dir, output_tag));
	fprintf(output, "TOTAL_NODES=%lu\n", my_graph.getNumNodes());
	fprintf(output, "TOTAL_EDGES=%lu\n", my_graph.getNumEdges());
	fclose(output);
	
	return 0;
}
