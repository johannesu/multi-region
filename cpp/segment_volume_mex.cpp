//
// Multi region segmentation using Lagrangian duality (and roof duality).
// Johannes Ulén and Petter Strandmark.
//
#define TIMING
#define PRINT_TIMING
#include "mex.h"
#include "maxflow-v3.02.src/graph.h"
#include "QPBO-v1.3.src/QPBO.h" 
#include "utils/mexutils.h"
#include "utils/mextiming.h"
#include "utils/cppmatrix.h"
#include <stdexcept>
#include <stdlib.h> 
#include <cmath>
#include <memory>
#include <algorithm>
#include <limits>

using std::min;

// Regions indices:
// Stacked by region i.e. 0--n-1, region 0, n--2*n-1 region 1, et c.
// Region 1 is background and no extra variables is needed hence r-1.
// i: variable number
// r: region number
// num_variables: number of variables
inline int get_index(int i, int r, int num_variables)
{
	return i +(r-1)*num_variables;
}

// n!
int factorial(int n)
{
	if (n <= 1)
		return 1;

	int fact = 1;

	for(int i=1; i <= n;i++)
			fact=fact*i;

	return(fact);
}

int nchoosek(int n, int k)
{
	return factorial(n) / (factorial(n-k)*factorial(k));
}

// Debug function
int count_infeasibilities(std::vector<signed char>& labelling,
				   				 std::vector< std::vector<int> >& exclusions,
									 int num_variables)
{
	int infeasible = 0;

	for (int i = 0; i< num_variables; i++)
	{
		for (int e = 0; e < exclusions.size(); e++)
		{
			int hits = 0;
			for (int c = 0; c < exclusions[e].size(); c++)
			{
				if (labelling[get_index(i,  exclusions[e][c], num_variables)] == 1)
					hits++;

				if (hits == 2)
				{
					infeasible++;
					break;
				}
			}
		}
	}	
	return infeasible;
}


void error_function(const char *msg) 
{
	throw runtime_error(msg);
}

void error_function(char *msg) 
{
	throw runtime_error(msg);
}

inline double reciprocal_pairwise_no_trunc(double pixel1, double pixel2, double beta, double tau)
{
	return   1 / ( 1 + beta * std::abs(pixel1 - pixel2));
}

inline double reciprocal_pairwise(double pixel1, double pixel2, double beta, double tau)
{
	 return  1 / ( 1 + beta * min(std::abs(pixel1 - pixel2), tau) ); 
}

inline double scalar_pairwise(double pixel1, double pixel2, double beta, double tau)
{
		return 1;
}

inline double power_pairwise_no_trunc(double pixel1, double pixel2, double beta, double tau)
{
	return std::pow( std::abs(pixel1 - pixel2), beta);
}


inline double power_pairwise(double pixel1, double pixel2, double beta, double tau)
{
	return std::pow(  min ( std::abs(pixel1 - pixel2), tau), beta);
}

double (*pairwise)(double pixel1, double pixel2, double beta, double tau) = scalar_pairwise;

// Valid ind
bool valid_index(int n0, int n1, int n2, int n3, std::vector<int>& dims) 
{
	if ( (n0 > dims[0]-1 || n1 > dims[1]-1 || n2 > dims[2]-1 || n3 > dims[3] -1) || (n0 < 0 || n1 < 0 || n2 < 0 || n3 < 0) ) 
		return false;
	else
		return true;
}

// Linear index
inline int linear_index(int n0, int n1, int n2, int n3, std::vector<int>& dims) 
{
	if (!valid_index(n0,n1,n2,n3, dims)) 
		mexWarnMsgTxt("Out of range subscript");
	else 
		return  n0 + n1*dims[0] + n2*dims[0]*dims[1] + n3*dims[0]*dims[1]*dims[2];
}

// Make a solution feasible.
void make_feasible(std::vector<signed char>& labelling,
				   				 std::vector< std::vector<int> >& exclusions,
									 std::vector<int>& inclusions,
									 int num_variables)
{	
	for (int e = 0; e < exclusions.size(); e++)
	{
		// Randomly keep one label when resolving exclusion conflicts.
		int keep  = (int) rand()%(exclusions[e].size());

		for (int i = 0; i< num_variables; i++)
		{	
			bool infeasible = false;

			// Making it possible pass unlabelled solutions (values with -1).
			int hit = 0;
			for (int c = 0; c < exclusions[e].size(); c++)
			{
				if (labelling[get_index(i,  exclusions[e][c], num_variables)] == 1)
					hit++;

				if (hit == 2)
				{
				 	infeasible = true;
				 	break;
				}
			}
			
			// Making feasible
			if (infeasible)
			{
				for (int c = 0; c < exclusions[e].size(); c++)
				{
					if (c == keep)
						continue;

					labelling[get_index(i,  exclusions[e][c], num_variables)] = 0;
				}
			}
		}
	}

	// It possible that the previous step have broken some inclusion constraint.
	// Goes through all inclusion constraints until they are all satisfied

	// Upper bound on needed loops is number of inclusions.
	for (int constraint = 0; constraint < inclusions.size(); constraint++ )
	{
		bool success = true;

		for (int i = 0; i< num_variables; i++)
		{
			for (int j = 1; j < inclusions.size(); j++)
			{

				// Background is always enforced
				if (inclusions[j] == 0)
					continue;

				int child =  get_index(i, j,						num_variables);
				int parent = get_index(i, inclusions[j],num_variables);

				// Force compliance
				if ( (labelling[parent] == 0) && (labelling[child] == 1) )
				{
					success = false;
					labelling[child] = 0;
				}
			}
		}

		// No broken inclusions
		if (success)
			break;
	}
}

// Energy of a feasible labelling
template<typename REAL>
REAL compute_primal_energy(std::vector<signed char>& labelling,
														 std::vector<matrix <REAL> >& unaries,
														 std::vector<int>& inclusions,
														 const matrix<double>& connectivity,
														 const matrix<double>& connectivity_weights,
														 std::vector<int>& dims,
														 double beta, double tau)
{
	int num_variables = unaries[0].numel();
	REAL energy = 0;

	// Unary costs
	for (int parent = 1; parent < unaries.size(); parent++)
	{
		int child = inclusions[parent];
		for (int i = 0; i < num_variables; i++)
		{
			int index = get_index(i, parent, num_variables);

			if (labelling[index] == 1)
				energy +=  unaries[parent](i) - unaries[child](i);
		}
	} 

	// Pairwise energy
	for (int i = 0; i < dims[0]; ++i) {
	for (int j = 0; j < dims[1]; ++j) {
	for (int k = 0; k < dims[2]; ++k) {
	for (int l = 0; l < dims[3]; ++l) {	
		for (int ind=0; ind < connectivity.N; ++ind) 
		{
			int i2,j2,k2,l2;

			i2 = i + connectivity(0,ind);
			j2 = j + connectivity(1,ind);

			if (connectivity.M > 2)
				k2 = k + connectivity(2,ind);
			else
				k2 = 0;

			if (connectivity.M > 3)
				l2 = l + connectivity(3,ind);
			else 
				l2 = 0;

			if (valid_index(i2,j2,k2,l2, dims))
			{
				int sourceid = linear_index(i,j,k,l, dims);
				int targetid = linear_index(i2,j2,k2,l2, dims);

				// Loop over regions
				for (int r = 1; r < unaries.size(); r++)
				{
					if ( (labelling[get_index(sourceid, r, num_variables)] == 0) && 
							 (labelling[get_index(targetid, r, num_variables)] == 1) )
					{
						REAL pixel1 = unaries[r](sourceid);
						REAL pixel2 = unaries[r](targetid);

						energy += connectivity_weights(ind)* pairwise(pixel1,pixel2, beta, tau);
					}
				}
			}
		}
	}
	}
	}
	}

	return energy;
}

template<typename REAL>
void segment_volume(int            nlhs, 		/* number of expected outputs */
								 mxArray        *plhs[],	/* mxArray output pointer array */
								 int            nrhs, 		/* number of inputs */
								 const mxArray  *prhs[]		/* mxArray input pointer array */)
{
	ASSERT(nlhs == 3);
	ASSERT(nrhs == 6);
	
	startTime();

	// Reading input
	ASSERT(mxIsCell(prhs[0]));
	ASSERT(mxIsCell(prhs[2]));

	// Background is counted
	int num_regions = mxGetNumberOfElements(prhs[0]);
	int num_exclusions = mxGetNumberOfElements(prhs[2]);

	ASSERT(num_regions > 0);

	vector<matrix <REAL> > unaries;
	vector< vector<int> > exclusions;
	
	for (int i = 0; i < num_regions; i++)
		 unaries.push_back(mxGetCell(prhs[0], i));

	std::vector<int> dims(4);

	dims[0] = unaries[0].M;
	dims[1] = unaries[0].N;
	dims[2] = unaries[0].O;
	dims[3] = unaries[0].P;
	
	const matrix<int> inclusions_matrix(prhs[1]);
	std::vector<int> inclusions(mxGetNumberOfElements(prhs[1]));
	
	// Check inclusions
	// Number of hard coded inclusions 
	int num_inclusions = 0;
	for (int i = 0; i < inclusions.size(); i++)
	{
		// Change index to 0 based index
		inclusions[i] = inclusions_matrix(i) -1;

		// Inclusion in background does not need any
		// extra edges
		if (inclusions[i] > 0)
			num_inclusions++;
	}

	//
	// Parse exclusion constraints
	for (int i = 0; i < num_exclusions; i++)
	{
		const matrix<int> exclusions_matrix(mxGetCell(prhs[2], i));
		std::vector<int> exclude_list(exclusions_matrix.numel());

		// Index starting at 0
		for (int j = 0; j < exclusions_matrix.numel(); j++)
			exclude_list[j] = exclusions_matrix(j) - 1;

		exclusions.push_back(exclude_list);
	}

	// To be used for roof duality
	// number of edges needed for the exclusion constraints for 1 pixel.
	int num_exclude_edges = 0;
	for (int e = 0; e < exclusions.size(); e++)
		num_exclude_edges += nchoosek(exclusions[e].size(),2);

	// ASSERT that the size of each unary cost is equal
	for (int i = 1; i < num_regions; i++)
	{
		ASSERT(unaries[0].M == unaries[i].M);
		ASSERT(unaries[0].N == unaries[i].N);
		ASSERT(unaries[0].O == unaries[i].O);
		ASSERT(unaries[0].P == unaries[i].P);
	}
	
	ASSERT(inclusions.size() == num_regions);

	// Number of voxels
	int num_variables =  dims[0]*dims[1]*dims[2]*dims[3];
	int num_binary_variables = num_variables*(num_regions-1);

	const matrix<double> connectivity(prhs[3]);
	const matrix<double> connectivity_weights(prhs[4]);
	int connectivity_size = connectivity.N;

	// One per variable.
	int num_include_edges = num_inclusions*num_variables;
	num_exclude_edges *= num_variables;

	ASSERT(connectivity.M > 1 && connectivity.M <= 4)
	ASSERT(connectivity_weights.numel() == connectivity_size);

	//Structure to hold and parse additional parameters
	MexParams params(1, prhs+(nrhs-1));

	// MATLAB defaults everything to double.
	double 	lambda = params.get<double>("lambda", 1.0);       
	double 	beta = params.get<double>("beta", 1.0);
	double  tau = params.get<double>("tau",-1.0);
	bool		verbose = params.get<bool>("verbose", true);  
	int			maxiter = params.get<double>("maxiter", 25); 
	REAL 	maxRelGap = params.get<double>("maxRelGap", 1e-4);
	bool 		langrian_duality 	= params.get<bool>("langrian_duality", false); 
	bool 		RD 	= params.get<bool>("rd", false);
	bool 		improve = params.get<bool>("improve", false);

	ASSERT(maxiter > 0);
	ASSERT(maxRelGap > 0);
	ASSERT(lambda > 0);

	// Incorporate the lambda into the weighting factor 
	for (int i = 0; i < connectivity_weights.numel(); i++)
			connectivity_weights.data[i] *= lambda;

	string str_stepsize_rule = params.get<string>("steprule","adaptive");
	enum Steprule {adaptive,oneoverk};
	Steprule steprule = adaptive;

	if (str_stepsize_rule == "adaptive") 
		steprule = adaptive;
	else if (str_stepsize_rule == "oneoverk") 
		steprule = oneoverk;
	else 
		throw runtime_error("Unknown step size rule");

	string str_reg_type = params.get<string>("regularization_type", "scalar");


	// Each regularization is defined in two version one with and one without truncation
	// this is slightly more efficient then setting tau =  std::numeric_limits<double>::infinity(); 
	if (tau < 0.0)
	{
		if (str_reg_type == "reciprocal")
			pairwise = &reciprocal_pairwise_no_trunc;
		else if (str_reg_type == "scalar")
			pairwise = &scalar_pairwise;
		else if (str_reg_type == "power")
			pairwise = &power_pairwise_no_trunc;
		else
			throw runtime_error("Unknown regularization type");
	} else
	{
		if (str_reg_type == "reciprocal")
			pairwise = &reciprocal_pairwise;
		else if (str_reg_type == "scalar")
			pairwise = &scalar_pairwise;
		else if (str_reg_type == "power")
			pairwise = &power_pairwise;
		else
			throw runtime_error("Unknown regularization type");
	}



	// Arbitrary label for unlabelled nodes
	int rd_default_label = 0;


	if (verbose)
	{
		endTime("# Settings: \n");

		if (RD)
		{
			endTime("Solver: Roof duality \n");
			if (improve)
					endTime("Using improve \n");
		} 
		else
		{	
			endTime("Solver: Lagrangian duality \n");
			endTime("Maximum number of iterations: %d \n", maxiter);
			endTime("Step size rule: %s \n", str_stepsize_rule.c_str());
		}

		endTime("Regularization type: %s \n", str_reg_type.c_str());
		endTime("Regularization lambda: %g \n", lambda);
		endTime("Regularization beta: %g \n", beta);

		if (tau < 0.0)
			endTime("Using no truncation (tau < 0.0).\n");
		else
			endTime("Regularization tau: %g  \n", tau);

		endTime("Verbose: %s \n", (verbose)?"true":"false");
		endTime("\n");
		endTime("# Input data: \n");;
		endTime("Image size: %d x %d x %d x %d \n", dims[0], dims[1], dims[2],dims[3]);
		endTime("Number of regions: %d \n", num_regions);
		endTime("Connectivity size: %d \n", connectivity_size);
		endTime("Number of inclusion constraints (excluding background inclusion): %d \n", num_inclusions);
		endTime("Number of exclusion constraints: %d \n", exclusions.size());
	}

	// Setting output
	matrix<signed char> labelling(num_binary_variables);
	matrix<REAL> largest_lower_bound(1);
	matrix<REAL> final_energy(1);
	vector< signed char > current_labelling(num_binary_variables,0);

	plhs[0] = labelling;
	plhs[1] = largest_lower_bound;
	plhs[2] = final_energy;

	for (int i = 0; i < num_binary_variables; i++)
		labelling(i) = 0;

	REAL inf = std::numeric_limits<REAL>::infinity();
	REAL min_feasible_energy = inf;
	REAL max_lower_bound = -inf;
	int unlabelled = 0;

	// Fast upper bound given by is gien by
	// num_connectivity_edges = connectivity_size*num_binary_variables;
	// to save memory we calculate the exact value
	int num_connectivity_edges = 0;
	for (int i = 0; i < dims[0]; ++i) {
	for (int j = 0; j < dims[1]; ++j) {
	for (int k = 0; k < dims[2]; ++k) {
	for (int l = 0; l < dims[3]; ++l) {	
		for (int ind=0; ind < connectivity.N; ++ind) 
		{
			int i2,j2,k2,l2;

			i2 = i + connectivity(0,ind);
			j2 = j + connectivity(1,ind);

			if (connectivity.M > 2)
				k2 = k + connectivity(2,ind);
			else
				k2 = 0;

			if (connectivity.M > 3)
				l2 = l + connectivity(3,ind);
			else 
				l2 = 0;

			if (valid_index(i2,j2,k2,l2, dims))
				num_connectivity_edges+= (num_regions-1);
		}
	}
	}
	}
	} 

	if (verbose)
		endTime("# Starting optimization");

	// Solve using Lagrangian duality
	if (!RD)
	{
		// Setup problem
		// num_nodes, upper bound on num edges + data term edges
		int max_num_edges = num_connectivity_edges + 
												num_include_edges;

		Graph<REAL, REAL, REAL> graph(num_binary_variables, max_num_edges, error_function);
		graph.add_node(num_binary_variables);

		//
		// Unary costs
		//
		for (int parent = 1; parent < unaries.size(); parent++)
		{
			int child = inclusions[parent];

			// Usual cost
			for (int i = 0; i < num_variables; i++)
			{
				REAL cost = unaries[parent](i) - unaries[child](i);
				graph.add_tweights(get_index(i, parent, num_variables), cost, 0);
			}
		} 

		int num_pairwise_added = 0;
		//
		// Inclusions
		//
		for (int parent = 1; parent < inclusions.size(); parent++)
		{
			int child = inclusions[parent];

			// No edges needed to include background.
			if (child == 0)
				continue;

			// Submodular constraints
			for (int i = 0; i < num_variables; i++)
			{			
				graph.add_edge(get_index(i, child, num_variables), 
											 get_index(i, parent, num_variables), 
											 inf,0);
			}
		}

		int added_edges = 0;
		//
		// Add connectivity
		//
		for (int i = 0; i < dims[0]; ++i) {
		for (int j = 0; j < dims[1]; ++j) {
		for (int k = 0; k < dims[2]; ++k) {
		for (int l = 0; l < dims[3]; ++l) {	
			for (int ind=0; ind < connectivity.N; ++ind) 
			{
				int i2,j2,k2,l2;

				i2 = i + connectivity(0,ind);
				j2 = j + connectivity(1,ind);

				if (connectivity.M > 2)
					k2 = k + connectivity(2,ind);
				else
					k2 = 0;

				if (connectivity.M > 3)
					l2 = l + connectivity(3,ind);
				else 
					l2 = 0;

				if (valid_index(i2,j2,k2,l2, dims))
				{
					int sourceid = linear_index(i,j,k,l, dims);
					int targetid = linear_index(i2,j2,k2,l2, dims);

					// Loop over regions
					for (int r = 1; r < num_regions; r++)
					{
						REAL pixel1 = unaries[r](sourceid);
						REAL pixel2 = unaries[r](targetid);
						REAL w = connectivity_weights(ind)* pairwise(pixel1,pixel2, beta, tau);
 
						graph.add_edge(get_index(sourceid, r, num_variables), 
													 get_index(targetid, r, num_variables), 
													 w, 0);

						added_edges++;
					}
				}
			}
		}
		}
		}
		} 

		mexPrintf("added edges: %d \n", added_edges);

		//
		// Start solving
		//

		// Allocate memory
		// One multiplier per exclusion
		vector< vector <REAL> > lagrangian_multipliers;
		vector< vector <signed char> >  supergradients;

		// Used for adaptive-steprule
		vector< vector <REAL> > individual_step_size;
		vector< vector <signed char> >  prev_diff;
		
		for (int e = 0; e < num_exclusions; e++)
		{
			lagrangian_multipliers.push_back( vector<REAL>(num_variables,0) );
			individual_step_size.push_back( vector<REAL>(num_variables,1) );
			prev_diff.push_back( vector<signed char>(num_variables,0) );
		} 

		// First iteration there is no flow to reuse
		bool reuse = false; 
		for (int iter = 0; iter < maxiter; iter++)
		{
			REAL lower_bound = graph.maxflow(reuse);

			// next iteration reuse flows
			reuse = true; 
			
			// \lambda *-1 cost is not added to the graph 
			// so we remove it here.
			for (int e = 0; e < num_exclusions; e++) 
			{
				for (int i = 0; i < num_variables; ++i)  
				{
						lower_bound -= lagrangian_multipliers[e][i];
				}
			}

			if (max_lower_bound < lower_bound)
				max_lower_bound = lower_bound;
														
			// Go through every pixel and calculate supergradient
			int infeasibilities = 0;

			//  Exclusion constraint: x(r1)+x(r2)+....x(rn) \leq 1
			for (int i = 0; i < num_variables; ++i) 
			{

				// One multiplier per variable and constraint
				for (int e = 0; e < exclusions.size(); e++)
				{
					signed char supg = -1;

					for (int c = 0; c < exclusions[e].size(); c++)
							supg += graph.what_segment(get_index(i,  exclusions[e][c], num_variables));
	

					if (supg > 0) 
						infeasibilities++;

					// Will always result in 0 step
					if (supg == 0)
						continue;

					REAL step; 

					// Petter CVPR 2010 step sizes
					if (steprule == adaptive) 
					{
						// Half step size if we switched
						// from feasible to non-feasible or 
						// vice versa.
						if (prev_diff[e][i]*supg < 0)
							individual_step_size[e][i] /= 2;

						prev_diff[e][i] = supg;
						step = individual_step_size[e][i]*supg;
					}

					// Simple 1/k
					if (steprule == oneoverk)
						step = supg/(iter+1);

					// Update multipliers
					lagrangian_multipliers[e][i] += step;

					// Project onto feasible set
					if (lagrangian_multipliers[e][i] < 0)
					{
							step -= lagrangian_multipliers[e][i];
							lagrangian_multipliers[e][i] = 0;
					}
		
					// Add the cost lambda*x1 + lambda*x2+ ... and marking the nodes
					for (int c = 0; c < exclusions[e].size(); c++)
					{
						graph.add_tweights(get_index(i,  exclusions[e][c], num_variables),
															 step,
															 0);

						graph.mark_node(get_index(i,  exclusions[e][c], num_variables));
					}
				}
			}

			// The initial solutions is feasible: no need to perform iterations
			if (iter == 0 && infeasibilities == 0)
			{
				if (verbose)
					endTime("Energy: %g Lower bound : %g. Initial solution is feasible, we have global optimum directly. \n", min_feasible_energy, max_lower_bound);

				for (int i = 0; i < num_binary_variables; ++i)
					labelling(i) = graph.what_segment(i);

				min_feasible_energy = lower_bound;

				break;
			}

			for (int i = 0; i < num_binary_variables; ++i)
				current_labelling[i] = graph.what_segment(i);

		
			if (infeasibilities > 0)
				make_feasible(current_labelling, exclusions, inclusions, num_variables);

			// Get feasible energy, upper bound
			REAL feasible_energy =  compute_primal_energy<REAL>
															(current_labelling, unaries, inclusions, 
															connectivity, connectivity_weights,
															dims, beta, tau);

			// Is this the best primal solution so far?
			// If so save the solution
			if (feasible_energy < min_feasible_energy) 
			{
				min_feasible_energy = feasible_energy;

				for (int i = 0; i < num_binary_variables; ++i)
							labelling(i) = current_labelling[i];
			}

			REAL relative_gap;
			
			if (min_feasible_energy != 0)				
			 relative_gap = (min_feasible_energy - max_lower_bound) / std::abs(min_feasible_energy);
			else
				relative_gap = (min_feasible_energy - max_lower_bound) / 1e-10;

			if (verbose)
				endTime("Iter: %d, Energy: %g Lower bound %g Rel. gap: %g Infeasibilities: %d  \n", 
									iter, min_feasible_energy, max_lower_bound, relative_gap, infeasibilities);

			// Check convergence criteria.
			if (relative_gap < maxRelGap)
			{
					if (verbose)
					{
						endTime("Sufficiently small relative duality gap found.\n");
						endTime("Relative gap: %g  Needed relative gap: %g.   \n", relative_gap, maxRelGap);
					}

					break;
			} 

			if ((iter == maxiter-1) && verbose)
				endTime("Reached maximum number of iterations %d \n", iter+1);
		}	
	}  
	
	// Solve using Roof duality
	if (RD)
	{
		// Using std::numeric_limits<double>::infinty is not compatible with the QPBO code.
		inf = 1e10;
 
 		int max_num_edges = num_connectivity_edges + 
												num_include_edges + 
												num_exclude_edges;

		QPBO<REAL> qpbo(num_binary_variables, max_num_edges, error_function);
		qpbo.AddNode(num_binary_variables);

		for (int parent = 1; parent < unaries.size(); parent++)
		{
			int child = inclusions[parent];

			// Usual cost
			for (int i = 0; i < num_variables; i++)
			{
				REAL cost = unaries[parent](i) - unaries[child](i);
				int index = get_index(i, parent, num_variables);
				
				qpbo.AddUnaryTerm(index, 0, cost);
			}
		} 

		// Inclusions; submodular
		for (int parent = 1; parent < inclusions.size(); parent++)
		{
			int child = inclusions[parent];

			// No edges needed to include background.
			if (child == 0)
				continue;

			for (int i = 0; i < num_variables; i++)
			{
				int source_id = get_index(i, child	, num_variables);
				int target_id = get_index(i, parent	, num_variables);
				
				qpbo.AddPairwiseTerm(source_id, target_id, 0, inf,0, 0);
			}
		}

		// Exclusions; Non submodular
		// Need to add all combinations
		for (int i = 0; i < num_variables; i++)
		{
			for (int e = 0; e < exclusions.size(); e++)
			{
				for (int c1 = 0; c1 < exclusions[e].size()-1; c1++) {
				for (int c2 = c1+1; c2 < exclusions[e].size(); c2++) {
		
					int source_id = get_index(i, exclusions[e][c1], num_variables);
					int target_id = get_index(i, exclusions[e][c2], num_variables);

					qpbo.AddPairwiseTerm(source_id, target_id, 0, 0, 0, inf);
				}
				}
			}
		} 

		// Add connectivity
		for (int i = 0; i < dims[0]; ++i) {
		for (int j = 0; j < dims[1]; ++j) {
		for (int k = 0; k < dims[2]; ++k) {
		for (int l = 0; l < dims[3]; ++l) {	
			for (int ind=0; ind < connectivity.N; ++ind) 
			{
				int i2,j2,k2,l2;

				i2 = i + connectivity(0,ind);
				j2 = j + connectivity(1,ind);

				if (connectivity.M > 2)
					k2 = k + connectivity(2,ind);
				else
					k2 = 0;

				if (connectivity.M > 3)
					l2 = l + connectivity(3,ind);
				else 
					l2 = 0;

				if (valid_index(i2,j2,k2,l2, dims))
				{
					int sourceid = linear_index(i,j,k,l, dims);
					int targetid = linear_index(i2,j2,k2,l2, dims);

					// Loop over regions
					for (int r = 1; r < num_regions; r++)
					{
						REAL pixel1 = unaries[r](sourceid);
						REAL pixel2 = unaries[r](targetid);
						REAL w = connectivity_weights(ind)* (REAL)pairwise(pixel1,pixel2, beta, tau);

						qpbo.AddPairwiseTerm(get_index(sourceid, r, num_variables), 
																 get_index(targetid, r, num_variables), 
																 0, w, 0, 0);
					}
				}
			}
		}
		}
		}
		} 

		// Solve
		//
		qpbo.MergeParallelEdges();
		qpbo.Solve();       
		qpbo.ComputeWeakPersistencies();

		// First solution
		for (int i = 0; i < num_binary_variables; ++i) 
		{
			current_labelling[i] = qpbo.GetLabel(i);

			if (current_labelling[i] < 0) 
			{
				unlabelled++;
				current_labelling[i] = rd_default_label;
			}
		}

		max_lower_bound = qpbo.ComputeTwiceLowerBound()/2;
		REAL qpbo_energy = qpbo.ComputeTwiceEnergy()/2;

		if (unlabelled > 0)
			make_feasible(current_labelling, exclusions, inclusions, num_variables);

		min_feasible_energy =  compute_primal_energy(current_labelling, unaries, inclusions,
																								 connectivity, connectivity_weights,
																								 dims, beta, tau);
		
		// Save feasible solution
		for (int i =0; i < num_binary_variables; i++)
			labelling(i) = current_labelling[i];


		if (verbose)
			endTime("Iter: %d, Energy: %g Lower bound %g unlabelled: %d  \n", 
									0, min_feasible_energy, max_lower_bound, unlabelled);

		if (unlabelled == 0 && verbose)
			endTime("Global optimum found. \n");

		// Improve
		if (unlabelled > 0 && improve)
		{
			endTime("Starting improve for a maximum of %d iterations.\n", maxiter);

			srand(0);  // QPBO-code asks user to reinitialize the random number
								// for reproducibility reasons we simply chose 0.

			bool improve_result;
			REAL relative_gap;
			for (int iter = 0; iter < maxiter; iter++)
			{
					improve_result = qpbo.Improve();

					if (improve_result)
					{
						REAL lower_bound = qpbo.ComputeTwiceLowerBound()/2;

						if (max_lower_bound < lower_bound)
							max_lower_bound = lower_bound;


						// Read solution
						for (int i = 0; i < num_binary_variables; ++i) 
						{
							current_labelling[i] = qpbo.GetLabel(i);

							if (current_labelling[i] < 0) 
								current_labelling[i] = rd_default_label;
						}

						make_feasible(current_labelling, exclusions, inclusions, num_variables);

						REAL feasible_energy =  compute_primal_energy(current_labelling, unaries, inclusions,
																														connectivity, connectivity_weights,
																														dims, beta, tau);
						// Save best solution so far.
						if (feasible_energy < min_feasible_energy) 	
						{
							min_feasible_energy = feasible_energy;

							for (int i = 0; i < num_binary_variables; ++i)
										labelling(i) = current_labelling[i];
						} else
						{
							// Improve may not improve on _feasible_ energy 
							improve = false;
						}
						
						
						if (min_feasible_energy != 0)				
							 relative_gap = (min_feasible_energy - max_lower_bound) / std::abs(min_feasible_energy);
						else
								relative_gap = (min_feasible_energy - max_lower_bound) / 1e-10;

					}

				endTime("Iter: % Energy: %g Lower bound %g Rel. gap \n", iter+1, min_feasible_energy, max_lower_bound, relative_gap);

				// Check convergence criteria.
				if (relative_gap < maxRelGap)
				{
					if (verbose)
					{
						endTime("Sufficiently small relative duality gap found.\n");
						endTime("Relative gap: %g  Needed relative gap: %g.   \n", relative_gap, maxRelGap);
					}
						break;
				} 
			}
		}
	}

	largest_lower_bound(0) = max_lower_bound;
	final_energy(0)  = min_feasible_energy;

	endTime("Done. \n");
}

void mexFunction(int            nlhs, 		/* number of expected outputs */
								 mxArray        *plhs[],	/* mxArray output pointer array */
								 int            nrhs, 		/* number of inputs */
								 const mxArray  *prhs[]		/* mxArray input pointer array */)
{

	ASSERT(mxIsCell(prhs[0]));
	mxArray * unary_ptr = mxGetCell(prhs[0], 0);
	
	if (mxGetClassID(unary_ptr) == mxSINGLE_CLASS)
		 	segment_volume<float>(nlhs, plhs, nrhs, prhs);
	else if (mxGetClassID(unary_ptr) == mxDOUBLE_CLASS)
		segment_volume<double>(nlhs, plhs, nrhs, prhs);
	else
		 mexErrMsgTxt("Only double and single/float unary costs supported");

}