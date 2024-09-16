#include "Graph.h"
#include "Vertex.h"
#include "Edge.h"
#include "Timer.h"
#include "assert.h"

#include <ilcplex/ilocplex.h>
#include <ilconcert/ilomodel.h>
#include <vector>
#include <iostream>

const int SEED = 0;
const int NUM_THREADS = 2;
const int WORK_MEM = 12000;
const int TIME_LIMIT = 18000; // Timit limit in seconds for CPLEX
const double DESNITY_THRESHOLD = 0.5;
double TOL = 0.000001;

// void buildModel_S(IloModel &, IloNumVarArray &, Graph &, const std::size_t, const std::size_t, const double);
void buildModel_S(IloModel &,
									IloNumVarArray &,
									IloNumVarArray &,
									IloNumVarArray &,
									Graph &,
									const std::size_t,
									const std::size_t,
									const double);
std::size_t GetIndexFromVertClust(const std::size_t, const std::size_t, const std::size_t);
void ConvertToClusterNum(std::vector<std::vector<bool>> &, const std::size_t , std::vector<int> &);
void RecordResult(const std::string &, const std::string &, std::vector<int> &, const std::size_t);
void RecordObjective(const std::string &, const std::size_t, const double, const double);
void ReadClusterOutputFile(const std::string &, const std::size_t , std::vector<int> &);
std::size_t get_array_idx(const std::size_t i, const std::size_t j, const std::size_t num_cols);

int main(int argc, char **argv) {
	// Check inputs
	if (argc != 6) {
		printf("ERROR: Incorrect Inputs.\n");
		printf("Usage: %s <input.gml> <output_dir> <output_tag> <total_nodes> comp_num\n", argv[0]);
		exit(1);
	}
	
	Graph my_graph;
	std::vector<int> clust_result, incumb_result;
	std::string gml_file, output_dir, runtag;
	std::size_t total_nodes, total_edges, comp_num;
	Timer timer;
	double objValue;
	bool use_incumbent = false;
		
	// Read in user input to local variables
	gml_file = argv[1];
	output_dir = argv[2];
	runtag = argv[3];
	total_nodes = atoi(argv[4]);
	comp_num = atoi(argv[5]);
	
	// Read in GML/Graph file
	my_graph.readFile(gml_file);	
	clust_result.resize(my_graph.getNumNodes(), 0);

	std::string obj_file_name = output_dir + "S_ObjValue.csv";
	const double n_sq_n = (double)(my_graph.getNumNodes() * (my_graph.getNumNodes() - 1));
	double density = 2*my_graph.getNumEdges() / n_sq_n;
	printf("Density (comp %lu): %f\n", comp_num, density);	

	if ( density >= DESNITY_THRESHOLD) {		
		std::vector<int> single_clust(my_graph.getNumNodes(), 1);
		RecordResult(output_dir, runtag, single_clust, my_graph.getNumEdges());
		RecordObjective(obj_file_name, comp_num, 0, timer.elapsed_cpu_time());
		return 0;
	}

	const std::size_t K = static_cast<std::size_t>(floor(my_graph.getNumNodes()/2));
	const std::size_t S = my_graph.getNumNodes();
	double R = n_sq_n / (double)my_graph.getNumEdges();
	
	std::vector<std::vector<bool>> in_clust;
  std::vector<bool> in_clust_row(K+1, false);
  for(std::size_t i = 0; i < my_graph.getNumNodes(); ++i) {
    in_clust.push_back(in_clust_row);
	}

	printf("Starting S_MIP...");
  IloEnv env;
  try {
		IloCplex cplex(env);
		IloModel model(env);       

		const std::size_t num_x = (K+1) * my_graph.getNumNodes();
		const std::size_t num_b = K * my_graph.getNumNodes();
		const std::size_t num_ce = my_graph.getNumEdges();
		std::size_t num_vars = num_x + num_b + num_ce;
		fprintf(stderr, "Number of variables: %lu\n", num_vars)	;

    IloNumVarArray x(env, num_x, 0, 1, ILOINT);
		IloNumVarArray b(env, num_b, 0, 1, ILOINT);
		IloNumVarArray ce(env, num_ce, 0, 1, ILOINT);
		IloNumArray xCopy(env, num_x);	

		// Build model
    buildModel_S(model, x, b, ce, my_graph, K, S, R);
		cplex.extract(model);

		cplex.setParam(IloCplex::Param::RandomSeed, SEED);
		cplex.setParam(IloCplex::Param::Threads, NUM_THREADS);		
		cplex.setParam(IloCplex::Param::TimeLimit, TIME_LIMIT);
		//cplex.setParam(IloCplex::Param::WorkMem, WORK_MEM);
		//cplex.setParam(IloCplex::Param::MIP::Strategy::File, 2);
		//cplex.setParam(IloCplex::Param::Emphasis::Memory, 1);
		cplex.setOut(env.getNullStream());
		// cplex.setParam(IloCplex::SimDisplay, 0);
		// cplex.setParam(IloCplex::MIPDisplay, 1);

		timer.restart();
		cplex.solve();
		timer.stop();
		printf("complete\n");

		if (cplex.getStatus() == IloAlgorithm::Infeasible) {
			printf("Infeasible Solution\n");			
		} else if (cplex.getStatus() == IloAlgorithm::Optimal || cplex.getStatus() == IloAlgorithm::Feasible) {		
			objValue = cplex.getObjValue();
      cplex.getValues(x, xCopy);
			printf("MIP obj value: %f\n", objValue);
                  
			for(std::size_t i = 0; i < my_graph.getNumNodes(); ++i) {
				for(std::size_t j = 0; j < K+1; ++j) {
					if (xCopy[i*(K+1) + j] >= 1 - TOL)
            in_clust[i][j] = true;
				}
			}

			ConvertToClusterNum(in_clust, K, clust_result);
	  }
  } catch (IloException& e) {
		std::cerr << "Concert exception caught: " << e << std::endl;
	} catch (...) {
		std::cerr << "Unknown exception caught" << std::endl;
	}
	env.end();

	RecordResult(output_dir, runtag, clust_result, my_graph.getNumEdges());

	double S_obj = 1 - (objValue / n_sq_n);
	S_obj = S_obj*my_graph.getNumNodes()/total_nodes;
		
	RecordObjective(obj_file_name, comp_num, S_obj, timer.elapsed_cpu_time());

	printf("S: %f\n", S_obj);
	printf("Run Time: %f\n", timer.elapsed_cpu_time());
}

void buildModel_S(IloModel &model,
									IloNumVarArray &x,
									IloNumVarArray &b,
									IloNumVarArray &c,
									Graph &g,
									const std::size_t K,
									const std::size_t S, const double R) {
  IloEnv env = model.getEnv();
  
	std::size_t num_constraints = 0;
  // Require each node to be assign to exactly one cluster
	for(std::size_t i = 0; i < g.getNumNodes(); ++i) {
		IloExpr x_sum(env);
		for (std::size_t j = 0; j < K+1; ++j)
			x_sum += x[get_array_idx(i,j,K+1)];
		model.add(x_sum == 1); // Eq. 3
		x_sum.end();
		++num_constraints;
	}		
    
	// Find size of each cluster, nk (except for singleton cluster)
	for (std::size_t j = 0; j < K; ++j)	{
    IloExpr x_sum(env);
		IloExpr b_sum(env);
		
		for(std::size_t i = 0; i < g.getNumNodes(); ++i)
			x_sum += x[get_array_idx(i,j,K+1)];
		for(std::size_t i = 0; i < S; ++i)
			b_sum += b[get_array_idx(i,j,K)];
		
    model.add(x_sum == b_sum); // Eq 4
		++num_constraints;
    x_sum.end();
		b_sum.end();
	}
    
	// Set ce=1 if edge is intercluster or in singleton cluster
	Vertex* cur_vert;
	std::size_t e_ind = 0;
	for(std::size_t u = 0; u < g.getNumNodes(); ++u) {
		cur_vert = g.getPtrToVertex(u);
		
		for(std::size_t i = 0; i < cur_vert->getDegree(); ++i) {
			std::size_t v = cur_vert->getTargetId(i);
			if(u < v) {
				for(std::size_t j = 0; j < K; ++j) {
					model.add(x[get_array_idx(u,j,K+1)] - x[get_array_idx(v,j,K+1)] <= c[e_ind]); // Eq 5	
					++num_constraints;
				}				
				model.add(x[get_array_idx(u,K,K+1)] + x[get_array_idx(v,K,K+1)] <= static_cast<IloNum>(2)*c[e_ind]); // Eq 6
				++num_constraints;
					
				++e_ind;
			}
		}		
	}

	// Add objective function
	IloExpr obj(env);
	for(std::size_t i = 0; i < S; ++i) {
		for(std::size_t j = 0; j < K; ++j) {
			obj += 2 * static_cast<IloNum>(i) * b[get_array_idx(i,j,K)];
		}
	}

	for(std::size_t m = 0; m < g.getNumEdges(); ++m) {
		obj += R * c[m];
	}

	model.add(IloMinimize(env, obj));
	
	fprintf(stderr, "Number of constraints: %lu\n", num_constraints);
}

void ConvertToClusterNum(std::vector<std::vector<bool>> &in_clust,
												 const std::size_t k,
												 std::vector<int> &clust_num) {
	std::size_t cur_clust = 1;
	std::size_t node_count;
	// Loop through all non-singleton cluster
	for(std::size_t j = 0; j < k; ++j) {
    node_count = 0;
		// Check each row
		for(std::size_t i = 0 ; i < clust_num.size(); ++i) {
			if(in_clust[i][j]) {
        ++node_count;
			}
		}

		if(node_count == 1) {
			for(std::size_t i = 0 ; i < clust_num.size(); ++i) {
				if(in_clust[i][j]) {
					clust_num[i] = -1;
				}
			}
		} else if(node_count > 1) {
			for(std::size_t i = 0 ; i < clust_num.size(); ++i) {
				if(in_clust[i][j]) {
					clust_num[i] = cur_clust;
				}
			}
			++cur_clust;
		}
	}

	for(std::size_t i = 0 ; i < clust_num.size(); ++i) {
		if(in_clust[i][k] == 1) {
			clust_num[i] = -1;
		}
	}
}

void RecordResult(const std::string &output_dir,
									const std::string &run_tag,
									std::vector<int> &clust_assign,
									const std::size_t m) {
	FILE* output;
	std::string file_name;
	file_name = output_dir + run_tag + "_S.out";
	
	if((output = fopen(file_name.c_str(), "w")) == NULL) {
		printf("ERROR - Could not open *_MIP.out file(%s).\n", file_name.c_str());
		exit(1);
	}

	std::size_t num_clusts = 1;
	std::size_t num_singles = 0;
	for(std::size_t i = 0; i < clust_assign.size(); ++i) {
		if(clust_assign[i] == -1)
			++num_singles;
		else if(clust_assign[i] > num_clusts)
			num_clusts = clust_assign[i];
	}
	
	fprintf(output, "%lu nodes %lu clusters %lu edges\n", clust_assign.size(), num_clusts, m);
	for(std::size_t i = 0; i < clust_assign.size(); ++i)
		fprintf(output, "%i ", clust_assign[i]);

	fclose(output);
}

void RecordObjective(const std::string & file_name,
										 const std::size_t comp_num,
										 const double objective,
										 const double time) {
	FILE* output;
		
	if((output = fopen(file_name.c_str(), "a+")) == NULL) {
		printf("ERROR in Q_MIP.cpp- Could not open objective score file(%s).\n", file_name.c_str());
		exit(1);
	}

	fprintf(output, "%lu,%f,%f\n", comp_num, objective, time);	

	fclose(output);
}

void ReadClusterOutputFile(const std::string& out_file,
													 const std::size_t num_nodes,
													 std::vector<int>& clust_num) {
	clust_num.resize(num_nodes);

	FILE* input;
	char tmp_str[50];
	tmp_str[0] = '\0';
	
	printf("%s\n", out_file.c_str());
	if((input = fopen(out_file.c_str(), "r")) == NULL) {
		printf("ERROR - Cound not open cluster output file.\n");
		exit(1);
	}

	for(std::size_t n = 0; n < num_nodes; ++n) {
		if(fscanf(input, "%s", tmp_str) > 0)
			clust_num[n] = atoi(tmp_str);
	}

	fclose(input);
}

std::size_t get_array_idx(const std::size_t i,
													const std::size_t j,
													const std::size_t num_cols) {
	if (j >= num_cols) {
		fprintf(stderr, "ERROR - get_array_idx - Trying to access element out of bounds.\n");
		exit(1);
	}
	return (i * num_cols + j);
}