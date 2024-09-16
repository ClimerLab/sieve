#include "Graph.h"
#include "Vertex.h"

#include "assert.h"
#include <cstring>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <limits>

Graph::Graph() :	num_edges(0),
									total_num_nodes(0),
									total_num_edges(0),
									weighted(false) {}

Graph::~Graph() {}

std::size_t Graph::getNumNodes() const {
  return vertices.size();
}

std::size_t Graph::getNumEdges() const {
  return num_edges;
}

std::size_t Graph::getTotalNumNodes() const {
	return total_num_nodes;
}

std::size_t Graph::getTotalNumEdges() const {
	return total_num_edges;
}

Vertex* Graph::getPtrToVertex(const std::size_t index) {
  assert(index < vertices.size());
  return &(vertices[index]);
}

void Graph::setTotalNumNodes(const std::size_t _num_nodes) {
	total_num_nodes = _num_nodes;
}

void Graph::setTotalNumEdges(const std::size_t _num_edges) {
	total_num_edges = _num_edges;
}

void Graph::readFile(const std::string &file_name) {
  // Check if the file has a '.gml' extension
	if (file_name.find(".gml", file_name.size() - 4) != std::string::npos) {
		readGML(file_name); // Read in file as .gml
	}
	// Check if the file has a '.graph' extension
	else if (file_name.find(".graph", file_name.size() - 6) != std::string::npos) {
		readGraph(file_name); // Read in file as .graph		
	} else {
    std::cerr << "ERROR - Input file type not recongnized." << std::endl;
    exit(1);
	}	
}

void Graph::readGML(const std::string &file_name) {
  // Declare local variable
	FILE *input;
	char tmp_str[50];
	std::size_t min_node, max_node, i_node, tmp_id, tmp_source, tmp_target;
	bool id_found, node_found, bracket_found, continuous_nodes, edge_found, source_found, target_found, edge_added;
	float tmp_weight;
	tmp_str[0] = '\0';

  std::size_t num_nodes = 0;
  std::vector<std::size_t> node_id_conv;

	// Try to open file
  if((input = fopen(file_name.c_str(), "r")) == NULL) {
    std::cerr << "ERROR in Graph::ReadGML - Could not open file (" << file_name.c_str() << ")." << std::endl;
		exit(1);
	}

	// Read file until start of graph is reahed
	while (true) {
		// Read in next string
		if (fscanf(input, "%s", tmp_str) <= 0) {
      std::cerr << "ERROR in Graph::ReadGML - Beginning of graph was not found. No data." << std::endl;
			fclose(input);
			exit(1);
		}

		// Check if end-of-file was reached
		if (feof(input)) {
      std::cerr << "ERROR in Graph::ReadGML - Beginning of graph was not found. EOF reached." << std::endl;
			fclose(input);
			exit(1);
		}

		// Check if the start of the graph was reached
		if (strncmp(tmp_str, "graph", 5) == 0) {
			break; // Break from while loop
		}
	}

	// Read in graph to determine number of nodes
	while (true) {
		// Read in next string
		if (fscanf(input, "%s", tmp_str) <= 0) {
      std::cerr << "ERROR in Graph::ReadGML - Number of nodes was not found." << std::endl;
			fclose(input);
			exit(1);
		}

		// Check if end-of-file was reached
		if (feof(input)) {
			break; // Break while loop
		}

		// Check if the start of a node was found
		if (strncmp(tmp_str, "node", 4) == 0) {
			++num_nodes;
		}

		// Check if start of edges was found
		// This should always occur after all the nodes are declared
		if (strncmp(tmp_str, "edge", 4) == 0) {
			break; // Break while loop
		}
	}

	// Move stream location to the start of the file
	fseek(input, 0, SEEK_SET);

	// Allocate memory for verticies and node ID conversion array
	vertices.resize(num_nodes);
	node_id_conv.resize(num_nodes);

	// Initialize variables for next step
	min_node = std::numeric_limits<std::size_t>::max();
	max_node = std::numeric_limits<std::size_t>::min();
	i_node = 0;
	node_found = false;
	bracket_found = false;
	id_found = false;
	// Loop until all nodes are read
	while (i_node < num_nodes) {
		// Check if end-of-file was reached
		if (feof(input)) {
			std::cerr << "ERROR in Graph::ReadGML - Node structure not formatted correctly." << std::endl;
      fclose(input);
			exit(1);
		}

		// Read in next string
		if (fscanf(input, "%s", tmp_str) <= 0) {
			std::cerr << "ERROR in Graph::ReadGML - Node structure not formated correctly.\n" << std::endl;
			fclose(input);
			exit(1);
		}

		// Check if 'node' was found
		if (strncmp(tmp_str, "node", 4) == 0) {
			node_found = true;
		}

		// Check if '[' was found following a node
		if (node_found && (strncmp(tmp_str, "[", 1) == 0)) {
			bracket_found = true;
		}

		// Check if ID was found when expected
		if (bracket_found && (strncmp(tmp_str, "id", 2) == 0)) {
			id_found = true; // Set flag			

			// Read in next string
			if (fscanf(input, "%s", tmp_str) <= 0) {
        std::cerr << "ERROR in Graph::ReadGML - Beginning of graph was not found. No ID value." << std::endl;
				fclose(input);
				exit(1);
			}

			// Store  ID value
			tmp_id = atoi(tmp_str);
		}

		// Check if ']' was found when expected
		if (id_found && (strncmp(tmp_str, "]", 1) == 0)) {
			// Store the value in the conversion arrray
			node_id_conv[i_node] = tmp_id;

			// Increment the counter
			i_node++;

			// Check if current ID is larger than max
			if (tmp_id > max_node)
				max_node = tmp_id;

			if (tmp_id < min_node)
				min_node = tmp_id;

			// Reset flags
			node_found = false;
			bracket_found = false;
			id_found = false;
		}
	}

	// Check if node IDs form a continous array (i.e. all values from min_node ... max_node are in the graph)
	bool continous_nodes = (max_node - min_node + 1 == num_nodes);

	// Initialize variables for next step
	edge_found = false;
	source_found = false;
	target_found = false;
	edge_added = false;
	tmp_source = 0;
	tmp_target = 0;
	tmp_weight = 1;

	// Loop until all edges are read
	while (true) {
		// Check if end-of-file was reached
		if (feof(input))
			break; // Break loop

		// Read in next line of file
		fscanf(input, "%s", tmp_str);

		// Check if edge was found
		if (strncmp(tmp_str, "edge", 4) == 0)
			edge_found = true; // Set flag

		// Check if '[' was found (start of edge information)
		if (edge_found && (strncmp(tmp_str, "[", 1) == 0))
			bracket_found = true; // Set flag

		// Check if source was found
		if (bracket_found && (strncmp(tmp_str, "source", 6) == 0)) {
			// Set flag
			source_found = true; // Set flag

			// Read in next string
			if (fscanf(input, "%s", tmp_str) <= 0) {
        std::cerr << "ERROR in Graph::ReadGML - Source was not found." << std::endl;
				fclose(input);
				exit(1);
			}

			// Convert string to int and store in temp variables
			tmp_source = atoi(tmp_str);
		}

		// Check if target was found
		if (bracket_found && (strncmp(tmp_str, "target", 6) == 0)) {
			// Set flag
			target_found = true;

			// Read in next string
			if (fscanf(input, "%s", tmp_str) <= 0) {
        std::cerr << "ERROR in Graph::ReadGML - Target was not found." << std::endl;
				fclose(input);
				exit(1);
			}

			// Convert to integer and store in tmp value
			tmp_target = atoi(tmp_str);
		}

		// Check if weight was found
		if (bracket_found && (strncmp(tmp_str, "weight", 6) == 0)) {
			// Read in next string
			if (fscanf(input, "%s", tmp_str) <= 0) {
        std::cerr << "ERROR in Graph::ReadGML - Error when trying to find edge weight." << std::endl;
				fclose(input);
				exit(1);
			}

			// Convert to float and store in temp vairalbe
			tmp_weight = atof(tmp_str);
		}

		// Check if ']' was found (end of edge information)
		if (source_found &&
				target_found &&
				(strncmp(tmp_str, "]", 1) == 0)) {
			// Check if a continous set of nodes was found
			if (continous_nodes) {
				// Update values based on value of lowest node (this forces the node IDs to start from 0)
				tmp_source -= min_node;
				tmp_target -= min_node;
			} else {
				std::cerr << "ERROR in Graph:ReadGML - Nodes do not form a continous array." << std::endl;
        fclose(input);
				exit(1);
			}

			// Try to add edge to graph
			edge_added = this->addEdge(tmp_source, tmp_target, tmp_weight);
			
			// Check if edge was successfully added
			if (!edge_added) {
				std::cerr << "ERROR in Graph::ReadGML - Add Edge failed for source " << tmp_source + min_node <<
                 		 " to target " << tmp_target + min_node << "." << std::endl;
			}

			// Check if graph is directed
			if (!graph_directed) {
				// Try to add edge to graph
				edge_added = this->addEdge(tmp_target, tmp_source, tmp_weight);

				// Check if edge was successfully added
				if (!edge_added) {
					std::cerr << "ERROR in Graph::ReadGML - Add Edge failed for source " << tmp_target + min_node <<
                     	 " to target " << tmp_source + min_node << "." << std::endl;
				}
			}

			// Increment the number of edges in graph
			++num_edges;

			// Reset flags
			edge_found = false;
			bracket_found = false;
			source_found = false;
			target_found = false;
		}
	}

	// Close file
	fclose(input);

	total_num_nodes = vertices.size();
	total_num_edges = num_edges;	
}

void Graph::readGraph(const std::string &file_name) {
  // Declare local variables
	std::size_t num_spaces, source_id, target_id, num_nodes_expected, num_nodes_found, num_edges_expected, num_edges_found, edges_added;
	bool weighted = false, edge_added;
	double weight;
	std::string tmp_str;
	std::istringstream iss;
	std::ifstream input;

	// Open the file
	input.open(file_name.c_str());

	// Check if file opened 
	if (!input) {
    std::cerr << "ERROR in Graph::ReadGraph - Could not open file (" <<
    			       file_name.c_str() << ")." << std::endl;		
		exit(1);
	}

	// Read the first line from the file
	std::getline(input, tmp_str);

	// Count the number of delimiters in the first line. Assume the strings are delimited by ' '. Expect 1 or 2 values.
	num_spaces = std::count(tmp_str.begin(), tmp_str.end(), ' ');

	// Check if the number of strings found is an expected value
	if ((num_spaces < 1) || (num_spaces > 2)) {
        std::cerr << "ERROR in Graph::ReadGraph - Header row has " << num_spaces+1
            			<< "strings, but only 2 or 3 are expected." << std::endl;
		exit(1);
	}

	iss.str(tmp_str);
	iss >> num_nodes_expected;
	iss >> num_edges_expected;
	if (num_spaces == 2)
		iss >> weighted;
		
	// Initialize the number of nodes and egdes found
	num_nodes_found = 0;
	num_edges_found = 0;

	// Read the next line from the input file untill the end is reached
	while (std::getline(input, tmp_str))
		++num_nodes_found; // Increment number
	
	// Check if the expected number of nodes matches the found number of nodes
	if (num_nodes_found != num_nodes_expected) {
    std::cerr << "ERROR in Graph::ReadGraph - Expected " << num_nodes_expected << "nodes, but found " << num_nodes_found << std::endl;
		exit(1);
	}

	// Allocate memory based on number of nodes found
  vertices.reserve(num_nodes_found);
	Vertex new_vert;
	for(std::size_t i = 0; i < num_nodes_found; ++i)
		vertices.push_back(new_vert);

	// Reset to second line in file
	input.clear();
	input.seekg(0);

	// Read first line
	std::getline(input, tmp_str);

	if (!weighted) {
		// The source are based on the row in the file. Initialize to zero.
		source_id = 0;
		while (std::getline(input, tmp_str)) {
			if (tmp_str.length() > 0) {
				// Count the number of delimiters
				num_spaces = std::count(tmp_str.begin(), tmp_str.end(), ' ');
				// Check if delimiter was found at end of string
				if (tmp_str.compare(tmp_str.size() - 1, tmp_str.size(), " ") == 0)
					num_spaces--;
				// Clear then update the istringstream
				iss.clear();
				iss.str(tmp_str);

				// Reset the number of edges added for this current vertex
				edges_added = 0;
				// Loop untill all edges are added
				while (edges_added < num_spaces + 1) {
					iss >> target_id;
					target_id--; // Account for file being 1 based, but array being zero based
					weight = 1;
					edge_added = this->addEdge(source_id, target_id, weight);
					num_edges_found++;
					edges_added++;

					// Check if edge was added
					if (!edge_added) {
						std::cerr << "Error in Graph::ReadGraph - Add Edge failed for source " << source_id << " to target " << target_id << std::endl;
						exit(1);
					}
				}
			}

			// Go to next nose as the source
			source_id++;
		}
	} else {
		// The source are based on the row in the file. Initialize to zero.
		source_id = 0;
		while (std::getline(input, tmp_str)) {
			// Count the number of delimiters
			num_spaces = std::count(tmp_str.begin(), tmp_str.end(), ' ');
			// Check if delimiter was found at end of string
			if (tmp_str.compare(tmp_str.size() - 1, tmp_str.size(), " ") == 0)
				num_spaces--;

			// Clear then update the istringstream
			iss.clear();
			iss.str(tmp_str);

			// Reset the number of edges added for this current vertex
			edges_added = 0;
			// Loop untill all edges are added
			while (edges_added < num_spaces / 2 + 1) {
				iss >> target_id;
				target_id--; // Account for file being 1 based, but array being zero based
				iss >> weight;
				edge_added = this->addEdge(source_id, target_id, weight);				
				num_edges_found++;
				edges_added++;

				// Check if edge was added
				if (!edge_added) {
					std::cerr << "Error in Graph::ReadGraph - Add Edge failed for source " << source_id << " to target " << target_id << std::endl;
					exit(1);
				}
			}

			// Go to next nose as the source
			source_id++;
		}
	}


	// Account for the fact that the .graph file has edges to and from the same source, but only counts it as 1
	num_edges_found = num_edges_found / 2;
	num_edges = num_edges_found;

	// Close file
	input.close();

	// Check if the number of edges found backed the expected value
	if (num_edges_found != num_edges_expected) {
		std::cerr << "ERROR in Graph::ReadGraph - Expected " << num_edges_expected << "nodes, but found " << num_edges_found << std::endl;
		exit(1);
	}

	total_num_nodes = vertices.size();
	total_num_edges = num_edges;	
}

bool Graph::addEdge(const std::size_t source_id,
										const std::size_t target_id,
										const double weight) {
  assert(source_id < vertices.size());
  assert(target_id < vertices.size());

  return vertices[source_id].addEdge(target_id, weight);
}

bool Graph::removeEdge(const std::size_t source_id,
											 const std::size_t target_id) {
	assert(source_id < vertices.size());
	assert(target_id < vertices.size());

	return vertices[source_id].removeEdgeById(target_id);
}

std::size_t Graph::getDegree(const std::size_t id) const {
	assert(id < vertices.size());

	return vertices[id].getDegree();
}

std::vector<std::size_t> Graph::getMinNeighbors(const std::size_t source_id,
																								const std::size_t target_id) {
  std::size_t min_id = source_id; // Initialize to source id

  // Set min_id to 'target_id' if 'target_id' has lower degree
  if(vertices[target_id].getDegree() < vertices[target_id].getDegree())
    min_id = target_id;

  return vertices[min_id].getIds();
}

std::vector<std::size_t> Graph::getMinSepFromNeighbors(const std::size_t source_id,
																											 const std::size_t target_id) {
	std::size_t u, v;
	if(vertices[source_id].getDegree() < vertices[target_id].getDegree()) {
		u = source_id;
		v = target_id;
	} else {
		u = target_id;
		v = source_id;
	}

	bool uv_adj = false;
	if(isAdjacent(u,v)) {
		uv_adj = true;
		removeEdge(u,v);
		removeEdge(v,u);
	}

	std::vector<std::size_t> neighbors = vertices[u].getIds();
	// Check if v is a neighbor of u. Remove from vector if found
	for(std::size_t i = 0; i < neighbors.size(); ++i) {
		if(neighbors[i] == v)
			neighbors.erase(neighbors.begin() + i);
	}

	// Run BFS from v until a node in 'neighbors' is reached. Continue to run BFS until no more nodes can be reached
	// Declare local variables
	std::size_t start_index = 0, end_index = 0, node_id;
	std::vector<std::size_t> queue;
	std::vector<bool> visited, is_neighbor;
	Vertex* cur_vert;
	
	// Allocate memory
	visited.resize(getNumNodes(), false);
	is_neighbor.resize(getNumNodes(), false);
	queue.resize(getNumNodes());

	// Add the v to the queue
	queue[end_index++] = v;
	visited[v] = true;

	// Loop through all neighbors of u
	for(std::size_t i = 0; i < neighbors.size(); ++i)
		is_neighbor[neighbors[i]] = true;
	
	// Loop until queue is empty
	while (start_index < end_index) {
		// Get node ID and dequeue
		node_id = queue[start_index++];

		// Get pointer to current vertex
		cur_vert = getPtrToVertex(node_id);

		// Loop through all edges in adjacency list
		for(std::size_t i = 0; i < cur_vert->getDegree(); ++i) {
			// Check if target id is a neighbor of u
			if(is_neighbor[cur_vert->getTargetId(i)]) {
				visited[cur_vert->getTargetId(i)] = true;
			}
			// Check if target id is unvisited
			else if(!visited[cur_vert->getTargetId(i)]) {
				// Add to queue
				queue[end_index++] = cur_vert->getTargetId(i);
				// Mark as visited
				visited[cur_vert->getTargetId(i)] = true;			
			}
		}
	}
	
	std::vector<std::size_t> min_neigh;
	for(std::size_t i = 0; i < neighbors.size(); ++i) {
		if(visited[neighbors[i]])
			min_neigh.push_back(neighbors[i]);
	}

	if(uv_adj) {
		addEdge(u,v,1.0);
		addEdge(v,u,1.0);
	}

	return min_neigh;
}

bool Graph::isAdjacent(const std::size_t i, const std::size_t j) const {
	assert(i < getNumNodes());
	assert(j < getNumNodes());

	if(vertices[i].isAdjacent(j) || vertices[j].isAdjacent(i))
		return true;
	else
		return false;
}

std::size_t Graph::seperateGraphIntoComponents(const std::string &out_dir,
																							 const std::string &output_tag) {
	Vertex *cur_vert;
	std::size_t start_index, end_index, node_id, cur_comp=0, largest_comp=0;
	std::size_t num1=0,num2=0,num3=0,num4_10=0,num11_100=0,num101_1000=0,num1000p=0;
	FILE *gml_file, *nn_file, *clust_file, *out_file;
	std::string nn_filename, gml_filename, clust_filename, out_filename;
	std::ostringstream ss;

	std::vector<bool> visited(getNumNodes(), false);
	std::vector<std::size_t> queue(getNumNodes());
	std::vector<std::size_t> node_id_converted(getNumNodes());
	std::vector<std::size_t> singletons(getNumNodes());

	// Open file for singletons and doubletons
	clust_filename = out_dir + output_tag + "_comp0.clust";
	clust_file = fopen(clust_filename.c_str(), "w");
	// Check that the file opened correctly
	if (clust_file == NULL) {
		printf("ERROR in Graph::seperateGraphIntoComponents - Could not open file: %s\n", clust_filename.c_str());
		exit(1);
	} else {
		fprintf(clust_file, "Cluster # (# nodes)\t|\tNodes\n");
	}

	// Loop through all nodes in the graph	
	for(std::size_t i = 0; i < getNumNodes(); ++i) {
		// Check if current node is a singleton
		if (getDegree(i) == 0) {
			fprintf(clust_file, "%lu (%i)\t\t\t\t|\t", num1,1);
			fprintf(clust_file, "%lu\n", i);
			singletons[num1++] = i;
			visited[i] = true;			
			if (largest_comp < 1)
				largest_comp = 1;
		}

		// Check if current node has been visisted
		else if(!visited[i]) {
			// Reset values
			start_index = 0;
			end_index = 0;

			// Set value in node ID conversion vector
			node_id_converted[i] = end_index;
			// Add noode to queue and increment end_index
			queue[end_index++] = i;
			// Set node as visited
			visited[i] = true;
			
			// Loop until queue is empty
			while(start_index < end_index) {
				// Get node ID and dequeue
				node_id = queue[start_index++];

				// Get pointer to current vertex
				cur_vert = getPtrToVertex(node_id);

				// Loop through all edges in adjacency list
				for(std::size_t i = 0; i < cur_vert->getDegree(); ++i) {
					// Check if target id is unvisited
					if(!visited[cur_vert->getTargetId(i)]) {
						// Set value in node ID conversion vector
						node_id_converted[cur_vert->getTargetId(i)] = end_index;
						// Add to queue
						queue[end_index++] = cur_vert->getTargetId(i);
						// Mark as visited
						visited[cur_vert->getTargetId(i)] = true;			
					}
				}
			}

			// Check if the current component is the largest
			if (end_index > largest_comp)
				largest_comp = end_index;

			// Increment current component number
			++cur_comp;
			ss.str(std::string()); // Clear stream contents
			ss << cur_comp; // Insert component number
				
			// Record *.nn file
			nn_filename = out_dir + output_tag + "_comp" + ss.str() + ".nn"; // Create full file name
			nn_file = fopen(nn_filename.c_str(),"w"); // Open file
			for (std::size_t j = 0; j < end_index; ++j)
				fprintf(nn_file, "%lu\n", queue[j]);			
			fclose(nn_file);

			// Record *.gml file
			gml_filename = out_dir + output_tag + "_comp" + ss.str() + ".gml";
			gml_file = fopen(gml_filename.c_str(), "w");

			fprintf(gml_file, "graph\n[");
			for (std::size_t j = 0; j < end_index; ++j)
				fprintf(gml_file, "\n  node\n  [\n    id %lu\n  ]", j);
			
			for (std::size_t j = 0; j < end_index; ++j) {
				// Get pointer to j-th vertex in queue
				cur_vert = getPtrToVertex(queue[j]);

				// Loop through all edges
				for(std::size_t k = 0; k < cur_vert->getDegree(); ++k) {
					std::size_t source_id = queue[j];
					std::size_t target_id = cur_vert->getTargetId(k);

					// Check if the converted source ID < converted target ID
					if(node_id_converted[source_id] < node_id_converted[target_id]) {
						if (this->weighted) {
							fprintf(gml_file, "\n  edge\n  [\n    source %lu\n    target %lu\n    weight %f\n  ]", node_id_converted[source_id], node_id_converted[target_id], cur_vert->getEdgeWeight(k));
						} else {
							fprintf(gml_file, "\n  edge\n  [\n    source %lu\n    target %lu\n  ]", node_id_converted[source_id], node_id_converted[target_id]);
						}
					}
				}
			}
			fprintf(gml_file, "\n]");
			fclose(gml_file);

			if(end_index == 2)
				++num2;
			else if (end_index == 3)
				++num3;
			else if (end_index < 11)
				++num4_10;			
			else if (end_index < 101)
				++num11_100;
			else if (end_index < 1001)
				++num101_1000;
			else
				++num1000p;			
		}
	}
	fclose(clust_file);


	// Wrtie nn file for comp 0
	std::string singleton_nn_filename = out_dir + output_tag + "_comp0.nn";
	FILE* singleton_nn_file= fopen(singleton_nn_filename.c_str(), "w");
	if(singleton_nn_file == NULL) {
		printf("ERROR - cound not open output file (%s)\n", singleton_nn_filename.c_str());
	}
	for(std::size_t i = 0; i < num1; ++i)
		fprintf(singleton_nn_file, "%lu\n", singletons[i]);
	fclose(singleton_nn_file);

	// Record results
	out_filename = out_dir + output_tag + ".stats";
	out_file = fopen(out_filename.c_str(), "w");
	fprintf(out_file, "Number of singletons: %lu\n", num1);
	fprintf(out_file, "Number of doubletons: %lu\n", num2);
	fprintf(out_file, "Number of tripletons: %lu\n", num3);
	fprintf(out_file, "Number of comps between 3-10: %lu\n", num4_10);
	fprintf(out_file, "Number of comps between 11-100: %lu\n", num11_100);
	fprintf(out_file, "Number of comps between 101-1000: %lu\n", num101_1000);
	fprintf(out_file, "Number of comps 1000+: %lu\n", num1000p);
	fprintf(out_file, "Largest component %lu nodes\n", largest_comp);
	fclose(out_file);

	return cur_comp;
}

void Graph::saveAsClusterAssignment(const std::string &output_filename) {
	Vertex *cur_vert;
	std::size_t start_index, end_index, node_id, cur_clust=0;
	FILE *out_file;

	std::vector<bool> visited(getNumNodes(), false);
	std::vector<std::size_t> queue(getNumNodes());
	std::vector<int> cluster_assingment(getNumNodes());
	std::vector<std::size_t> singletons(getNumNodes());

	// Loop through all nodes in the graph	
	for(std::size_t i = 0; i < getNumNodes(); ++i) {
		// Check if current node is a singleton
		if (getDegree(i) == 0) {
			cluster_assingment[i] = -1;
			visited[i] = true;
		}

		// Check if current node has been visisted
		else if(!visited[i]) {
			// Reset values
			start_index = 0;
			end_index = 0;

			cur_clust++;
			// Set cluster assignment number
			cluster_assingment[i] = cur_clust;
			// Add noode to queue and increment end_index
			queue[end_index++] = i;
			// Set node as visited
			visited[i] = true;
			
			// Loop until queue is empty
			while(start_index < end_index) {
				// Get node ID and dequeue
				node_id = queue[start_index++];

				// Get pointer to current vertex
				cur_vert = getPtrToVertex(node_id);

				// Loop through all edges in adjacency list
				for(std::size_t i = 0; i < cur_vert->getDegree(); ++i) {
					// Check if target id is unvisited
					if(!visited[cur_vert->getTargetId(i)]) {
						// Set value in node ID conversion vector
						cluster_assingment[cur_vert->getTargetId(i)] = cur_clust;
						// Add to queue
						queue[end_index++] = cur_vert->getTargetId(i);
						// Mark as visited
						visited[cur_vert->getTargetId(i)] = true;			
					}
				}
			}
		}
	}
	
	// Record results
	if((out_file = fopen(output_filename.c_str(), "w")) == NULL) {
		printf("ERROR in Graph::saveAsClusterAssignment - Could not open output file (%s)\n", output_filename.c_str());
		exit(1);
	}

	fprintf(out_file, "%lu nodes %lu clustes %lu edges\n", getNumNodes(), cur_clust, getNumEdges());

	for (std::size_t i = 0; i < getNumNodes(); ++i)
		fprintf(out_file, "%i ", cluster_assingment[i]);

	fclose(out_file);
}