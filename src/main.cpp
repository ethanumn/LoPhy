#include "utils.h"
#include "search.h"
#include "output.h"
#include "data_manager.h"
#include "csv.h"
#include "input.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <numeric>
#include <map>

using std::ifstream;
using std::getline;
using std::istringstream;
using std::stringstream;
using std::numeric_limits;
using std::map;
using std::stoi;

void read_parameters(int argc, char* argv[], Params & params);
void check_character_matrix(const vector<vector<int>> & character_matrix);
vector<string> get_cell_names(const string & cell_file, const int & num_cells);
void check_params(Params & params);
template<typename T>
void verify_matrix(const vector<vector<T>> & matrix, size_t num_rows=0, size_t num_columns=0);
vector<int> read_cell_samples(const string & cell_samples_file, const int & num_cells);

int main(int argc, char* argv[])
{

    // start time
    clock_t begin=clock();

    Params params;

	// read parameters
	read_parameters(argc, argv, params);

    // make sure we have everything we need to run scOrchard
    check_params(params);

    // read in character matrix 
    DataMatrix<int> character_data;
    read_csv(params.character_matrix_file, character_data, 0, 0);
    vector<vector<int>> character_matrix = character_data.get_data();
    const auto cell_names = character_data.get_index();
    check_character_matrix(character_matrix); // check character matrix to make sure it's valid
    verify_matrix(character_matrix);
    const size_t num_cells = character_matrix.size();
    const size_t num_variants = character_matrix[0].size();

    // read in variant reads
    DataMatrix<int> read_data;
    read_csv(params.variant_reads_file, read_data, 0, 0);
    vector<vector<int>> variant_reads = read_data.get_data();
    verify_matrix(variant_reads, num_cells, num_variants);

    // read in total reads
    read_data.clear();
    read_csv(params.total_reads_file, read_data, 0, 0);
    vector<vector<int>> total_reads = read_data.get_data();
    verify_matrix(total_reads, num_cells, num_variants);

    // read in region reads
    read_data.clear();
    read_csv(params.region_reads_file, read_data, 0, 0);
    vector<vector<int>> region_reads = read_data.get_data();
    vector<int> chromosomes = extract_chromosomes(read_data.get_columns());
    vector<string> region_names = extract_region_names(read_data.get_columns());

    // extract meta data
    DataMatrix<string> meta_data;
    read_csv(params.meta_file, meta_data, 0, 0);
    vector<int> locus_regions = convert_to_uint32(meta_data.get_column("REGION_INDEX"));
    vector<bool> is_snp = convert_to_bool(meta_data.get_column("SNP"));
    vector<bool> is_snv = convert_to_bool(meta_data.get_column("SBS"));
    vector<int> locus_samples = convert_to_uint32(meta_data.get_column("SAMPLE"));
    vector<string> gene_names = meta_data.get_column("NAME");
    gene_names.insert(gene_names.begin(), "root");

    // read in cell samples
    vector<int> cell_samples = read_cell_samples(params.cell_samples_file, num_cells);
    // const auto cell_names = get_cell_names(params.cell_file, num_cells);
    const auto path_with_prefix = out_path_prefix(params.character_matrix_file, params.out_path, params.out_prefix);

    // compute loss groups, only happens if K > 0 and there's a copy profile clustering
    data_manager manager(params.fp, 
                         params.fn, 
                         params.fn_indel,
                         params.hom_precision, 
                         params.het_precision, 
                         params.hom_precision_indel, 
                         params.het_precision_indel, 
                         params.dropout_concentration, 
                         params.dropout_rate_prior,
                         params.theta,
                         params.iters,
                         params.delta,
                         params.CNA_penalty,
                         params.tau1,
                         params.tau2);
    if(variant_reads.size() > 0 && total_reads.size() > 0)
        manager.fill(character_matrix,
                     variant_reads, 
                     total_reads, 
                     region_reads, 
                     cell_samples, 
                     locus_regions, 
                     is_snp, 
                     is_snv,
                     locus_samples, 
                     chromosomes);
    else
        throw invalid_argument("variant or total read counts are empty!");

    // initialize cache for scoring
    SCORE_CACHE cache;

    /* 
    Step 1
    -------
    Estimate sample-specific weights for the fraction of reads that fall into each region
    This is done by running the algorithm on each sample separately without considering CNAs
    */
    Tree tree;
    const auto & samples = manager.get_samples();
    for(const auto & sample : samples)
    {
        cout << "Step 1: build SNV only tree for sample " << sample << endl;
        // initialize Q with an empty tree
        tree = Tree(manager.get_num_loci(), num_cells, manager.get_dropout_rate_prior());

        // perform searchs
        tree = search(tree,
                      manager,
                      cache,
                      sample,
                      false,
                      true,
                      params.seed);

        manager.compute_region_weights(sample, tree.get_cell_assignments(), tree.get_parents());
        cout << endl << endl;

    }
    cout << endl << endl; 
    /*
    Step 2
    -------
    Construct tree containing SNVs, SNPs, and CNAs while adhering to longitudinal constraints
    */
    tree = Tree(manager.get_num_loci(), num_cells, manager.get_dropout_rate_prior());
    for(const auto & sample : samples)
    {
        cout << "Step 2: growing tree with data from sample " << sample << endl;
        cache.clear();
        tree = search(tree,
                      manager,
                      cache,
                      sample,
                      true,
                      false,
                      params.seed);
        tree.reset_dropout_rates(manager.get_dropout_rate_prior());
        cout << endl << endl;
    }

    // end time
    clock_t end=clock();
    double diffticks=end-begin;
    double diffms=(diffticks*1000)/CLOCKS_PER_SEC;

    cout << "Time elapsed: " << diffms << " ms"<< endl;
    write_all_results(path_with_prefix, 
                      cell_names,
                      gene_names,
                      region_names,
                      tree, 
                      manager,
		      params);

}

void print_args(void) {
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "--------------------------------------------  LoPhy input arguments  -------------------------------------------------------------------" << endl;
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << endl;
    cout << "Example call:                                                                                                                               " << endl;
    cout << "./LoPhy -c character_matrix.csv -r region_reads.csv -v variant_reads.csv -t total_reads.csv -m meta.csv -s cell_samples.txt -o outputs                                                          " << endl;
    cout << endl;
    cout << "------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "   Flag            |   Required   |  Description                                                                                            " << endl;
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    cout << "-c                 |     Yes      |  Cell x Mutation binary genotype matrix                                                                 " << endl;
    cout << "-o                 |     Yes      |  Path to write outputs to                                                                               " << endl;
    cout << "-m                 |     Yes      |  Meta data file containing information about each mutation                                              " << endl;
    cout << "-v                 |     Yes      |  Cell x Mutation variant read count matrix                                                              " << endl;
    cout << "-t                 |     Yes      |  Cell x Mutation total read count matrix                                                                " << endl;
    cout << "-r                 |     Yes      |  Region x Cell read depth matrix                                                                        " << endl;
    cout << "-s                 |     Yes      |  File containing which sample each cell is originally from                                              " << endl;
    cout << "-p                 |     No       |  Prefix to place on all output file                                                                     " << endl;
    cout << "-fp                |     No       |  Expected false positive rate (Default = 0.02)                                                          " << endl;
    cout << "-fn                |     No       |  Expected false negative rate for single base substitutions (Default = 0.02)                            " << endl;
    cout << "-fn-indel          |     No       |  Expected false negative rate for insertions and deletions (Default = 0.02)                             " << endl;
    cout << "-homp              |     No       |  Beta-binomial precision for SNVs when homozygous (Default = 50.0)                                      " << endl;
    cout << "-hetp              |     No       |  Beta-binomial precision for SNVs when heterozygous (Default = 8.0)                                     " << endl;
    cout << "--homp-indel       |     No       |  Beta-binomial precision for SNPs when homozygous (Default = 15.0)                                      " << endl;
    cout << "-hetp-indel        |     No       |  Beta-binomial precision for SNPs when heterozygous (Default = 4.0)                                     " << endl;
    cout << "-dropoutc          |     No       |  Concentration parameter for dropout likelihood (Default = 100.0)                                       " << endl;
    cout << "-dropoutp          |     No       |  Dropout rate prior probability (Default = 0.05)                                                        " << endl;
    cout << "-theta             |     No       |  Inverse-dispersion parameter for Negative binomial (Default = 6.0)                                     " << endl;
    cout << "-iters             |     No       |  Number of hill climbing iterations to perform for each subtree (Default = 200)                         " << endl;
    cout << "-delta             |     No       |  Maximum number of hill climbing iterations to perform without likelihood improvement (Default = 20)    " << endl;
    cout << "-tau1              |     No       |  Probability of performing an SNV relocation (Default = 0.5)                                            " << endl;   
    cout << "-tau2              |     No       |  Probability of adding a CNA clone (Default = 0.1)                                                      " << endl;   
    cout << "-penalty           |     No       |  Penalty for adding a CNA (Default = 85.0)                                                              " << endl;   
    cout << "-seed              |     No       |  Random seed (Default = random)                                                                         " << endl;
    cout << "--------------------------------------------------------------------------------------------------------------------------------------------" << endl;
    exit(0);
}

void check_params(Params & params) {
    bool passed_check = true;
    if(!std::filesystem::exists(params.character_matrix_file)) {
        cout << "Error with argument -c " << params.character_matrix_file << ". Path does not exist!" << endl;
        passed_check = false;
    } 
    if(!std::filesystem::exists(params.out_path)) {
        std::filesystem::create_directory(params.out_path);
        cout << "Created path " << params.out_path << endl;
    }
    if(!std::filesystem::exists(params.variant_reads_file)) {
        cout << "Error with argument -v " << params.variant_reads_file << ". Path does not exist!" << endl;
        passed_check = false;
    } 
    if(!std::filesystem::exists(params.total_reads_file)) {
        cout << "Error with argument -t " << params.total_reads_file << ". Path does not exist!" << endl;
        passed_check = false;
    } 
    if(!std::filesystem::exists(params.meta_file)) {
        cout << "Error with argument -m " << params.meta_file << ". Path does not exist!" << endl;
        passed_check = false;
    } 
    if(!std::filesystem::exists(params.region_reads_file)) {
        cout << "Error with argument -r " << params.region_reads_file << ". Path does not exist!" << endl;
        passed_check = false;
    } 
    if(params.fn < 0 || params.fn > 1.0) {
        cout << "Error with argument -fn " << params.fn << ". Must be a float in the range (0,1)." << endl;
        passed_check = false;
    } 
    if(params.fp < 0.0 || params.fp > 1.0) {
        cout << "Error with argument -fp " << params.fp << ". Must be a float in the range (0,1)." << endl;
        passed_check = false;
    } 
    if(params.tau1 < 0.0 || params.tau1 > 1.0) {
        cout << "Error with argument -tau1 " << params.fp << ". Must be a float in the range (0,1)." << endl;
        passed_check = false;
    } 
    if(params.tau2 < 0.0 || params.tau2 > 1.0) {
        cout << "Error with argument -tau2 " << params.fp << ". Must be a float in the range (0,1)." << endl;
        passed_check = false;
    } 
    if(params.iters < 0) {
        cout << "Error with argument -iters " << params.iters << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(params.delta < 0) {
        cout << "Error with argument -delta " << params.delta << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(params.hom_precision < 0) {
        cout << "Error with argument -homp " << params.hom_precision << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(params.het_precision < 0) {
        cout << "Error with argument -hetp " << params.het_precision << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(params.hom_precision_indel < 0) {
        cout << "Error with argument -homp-indel " << params.hom_precision_indel << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(params.het_precision_indel < 0) {
        cout << "Error with argument -hetp-indel " << params.het_precision_indel << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(params.dropout_rate_prior < 0.0 || params.dropout_rate_prior > 1.0) {
        cout << "Error with argument -dropoutrp " << params.dropout_rate_prior << ". Must be a float in the range (0,1)." << endl;
        passed_check = false;
    } 
    if(params.theta < 0.0) {
        cout << "Error with argument -theta " << params.theta << ". Must be non-negative." << endl;
        passed_check = false;
    } 
    if(params.CNA_penalty < 0) {
        cout << "Error with argument -penalty " << params.CNA_penalty << ". Must be a non-negative integer." << endl;
        passed_check = false;
    } 
    if(!passed_check)
    {
        cout << "Aborting scOrchard!" << endl;
        exit(1);
    }
}
template<typename T>
void verify_matrix(const vector<vector<T>> & matrix, size_t num_rows, size_t num_columns)
{
    if(matrix.empty())
        throw invalid_argument("One or more input matrices are empty!");

    if(num_rows > 0 && num_columns > 0)
    {
        if((matrix.size() != num_rows) || (matrix[0].size() != num_columns))
            throw invalid_argument("Input matrices sizes do not match!");

    }
}

// Helper function to find the index of an element in a list
int find_index(const std::vector<std::string>& reference_list, const std::string& element) {
    auto it = std::find(reference_list.begin(), reference_list.end(), element);
    if (it != reference_list.end()) {
        return std::distance(reference_list.begin(), it); // Return the index
    } else {
        return -1; // Element not found
    }
}

void check_character_matrix(const vector<vector<int>> & character_matrix)
{
    // check and make sure character_matrix only contains 0,1,3 as values 
    for(const auto & row : character_matrix)
        for(const auto val : row)
        {
            if(val == 0 || val == 1 || val == -1)
                continue;
            else
                goto failure_state;
        }
        
    // if completed return
    return; 

    // if failed, print error and exit
    failure_state: 
        cout << "character matrix input contains bad values! Each entry must be {-1, 0, 1}." << endl; 
        exit(1);
}

vector<int> read_cell_samples(const string & cell_samples_file, 
                                   const int & num_cells)
{
	vector<int> cell_samples;
	ifstream file(cell_samples_file.c_str());

    // if file doesn't exist then we'll assume all cells are in the same sample
	if (!file) 
    {
	    for(auto i = 0; i < num_cells; i++)
	    	cell_samples.push_back(0);
	} 
    else 
    {
        // extract each of the cell names from the file
        string samp;
        for (auto i = 0; i < num_cells; i++) {
            file >> samp;
            cell_samples.push_back(stoi(samp));
        }
    }

    // check to make sure they're contiguous, otherwise throw an error
    unordered_set<int> seen;
    int prev_sample = -1;

    for (int sample : cell_samples) {
        if (sample != prev_sample) {
            if (seen.count(sample)) {
                throw std::runtime_error("Data for each sample must appear contiguously in the input files.");
            }
            seen.insert(sample);
            prev_sample = sample;
        }
    }

    return cell_samples;
}

/* Reads cell names from a file */
vector<string> get_cell_names(const string & cell_file, 
                              const int & num_cells)
{

	vector<string> cell_names;
	ifstream file(cell_file.c_str());

    // if file doesn't exist then we'll just use integers for the gene names
	if (!file) {
	    for(auto i = 0; i < num_cells; i++)
	    	cell_names.push_back("s" + to_string(i));
	} 
    else 
    {
        // extract each of the cell names from the file
        string cell;
        for (auto i = 0; i < num_cells; i++) {
            file >> cell;
            cell_names.push_back(cell);
        }
    }

    return cell_names;
}

void read_parameters(int argc, char* argv[], Params & params){

	for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
            print_args();
        } else if (strcmp(argv[i], "-c") == 0) {
			if (i + 1 < argc) {params.character_matrix_file = argv[++i];}
        } else if (strcmp(argv[i], "-m") == 0) {
			if (i + 1 < argc) {params.meta_file = argv[++i];}
		} else if(strcmp(argv[i], "-o")==0) {
			if (i + 1 < argc) {params.out_path = argv[++i];}
		} else if(strcmp(argv[i], "-p")==0) {
			if (i + 1 < argc) {params.out_prefix = argv[++i];}
        } else if (strcmp(argv[i], "-v") == 0) {
			if (i + 1 < argc) {params.variant_reads_file = argv[++i];}
        } else if (strcmp(argv[i], "-t") == 0) {
			if (i + 1 < argc) {params.total_reads_file = argv[++i];}
        } else if (strcmp(argv[i], "-r") == 0) {
			if (i + 1 < argc) {params.region_reads_file = argv[++i];}
        } else if (strcmp(argv[i], "-s") == 0) {
			if (i + 1 < argc) {params.cell_samples_file = argv[++i];}
		} else if(strcmp(argv[i], "-fp")==0) {
			if (i + 1 < argc) {params.fp = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-fn")==0) {
			if (i + 1 < argc) {params.fn = atof(argv[++i]);}
		} else if(strcmp(argv[i],"-fn-indel")==0) {
			if (i + 1 < argc) {params.fn_indel = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-homp")==0) {
			if (i + 1 < argc) {params.hom_precision = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-hetp")==0) {
			if (i + 1 < argc) {params.het_precision = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-homp-indel")==0) {
			if (i + 1 < argc) {params.hom_precision_indel = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-hetp-indel")==0) {
			if (i + 1 < argc) {params.het_precision_indel = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-dropoutc")==0) {
			if (i + 1 < argc) {params.dropout_concentration = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-dropoutp")==0) {
			if (i + 1 < argc) {params.dropout_rate_prior = atof(argv[++i]);}
        } else if (strcmp(argv[i], "-theta") == 0) {
			if (i + 1 < argc) {params.theta = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-iters")==0) {
			if (i + 1 < argc) {params.iters = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-delta")==0) {
			if (i + 1 < argc) {params.delta = atoi(argv[++i]);}
		} else if(strcmp(argv[i], "-tau1")==0) {
			if (i + 1 < argc) {params.tau1 = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-tau2")==0) {
			if (i + 1 < argc) {params.tau2 = atof(argv[++i]);}
		} else if(strcmp(argv[i], "-penalty")==0) {
			if (i + 1 < argc) {params.CNA_penalty = atof(argv[++i]);}
		}  else if(strcmp(argv[i],"-seed")==0) {
			if (i + 1 < argc) {params.seed = atoi(argv[++i]) % numeric_limits<int>::max();}
		} 
	}
}
