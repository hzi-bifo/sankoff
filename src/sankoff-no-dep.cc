#include "tree.h"
#include "GetPot.cpp"

template<typename Value>
class blocked_range {
public:
    using const_iterator = Value;

    using size_type = std::size_t;

    blocked_range( Value begin_, Value end_) :
        my_end(end_), my_begin(begin_) {
    }

    const_iterator begin() const { return my_begin; }

    const_iterator end() const { return my_end; }

    size_type size() const {
        return size_type(my_end-my_begin);
    }

    //! True if range is empty.
    bool empty() const { return !(my_begin<my_end); }


private:
    Value my_end;
    Value my_begin;

};

template<typename V, typename F>
void parallel_for(const blocked_range<V>& r, const F& f) {
	f(r);
}

#include "sankoff.h"

int main(int argc, char* argv[]) {

	GetPot   vm(argc, argv);

	
// 	po::options_description desc("Fast Sankoff algorithm with parallel options.");
// 	desc.add_options()
// 	    ("help", "print help")
// 	    ("tree", po::value<string>(), "input tree file.")
// 	    ("nexus", po::value<string>(), "input nexus file.")
// 	    ("aln", po::value<string>()->required(), "aligned sequence file")
// 	    ("cost", po::value<string>(), "cost function file name")
// 	    ("cost-identity-aa", po::value<vector<int>>(), "Cost matrix is for chars of amino acids in addition to X and -. The parameters are cost of A-B where A and B are either amino acids (aa) X or gap with following order [identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap].")
// 	    ("cost-identity-dna", po::value<vector<int>>(), "Cost matrix is for chars of nucleic acid in addition to X and -. The parameters are cost of A-B where A and B are either nucleid acids (na) X or gap with following order [identical na] [non-identical na] [na-X] [X-X] [gap-X] [na-gap] [gap-gap].")
// 		//TODO: some default and simple cost functions
// //	    ("print-alignment", po::value<string>(), "print alignment")
// 	    ("omit-leaf-mutations", "omit mutations happen at leaf nodes")
// 	    ("nthread", po::value<int>()->default_value(oneapi::tbb::this_task_arena::max_concurrency()), "change number of default threads")
// 	    ("induce_tree_over_samples", po::value<bool>()->default_value(true), "remove nodes not in the alignment")
// 	;

	// po::variables_map vm;
	// po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.search("--help"))
	{
		cout << "Usage: " << argv[0] << " [options] ...\n";
		// cout << desc;
		//TODO: usage here
		return 0;
	}

	// oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, vm["nthread"].as<int>());
	// tick_count t0 = tick_count::now();

	//VM vm(argv);

	//Build tree from MSA removing excess nodes
	INode phylo;
	if (vm.search("--tree")) {
		ifstream fi(vm.next(""));
		phylo = load_tree<INode>(fi);
	} else 	if (vm.search("--nexus")) {
		ifstream fi(vm.next(""));
		phylo = load_nexus_tree<INode>(fi);
	} else {
		throw runtime_error("No tree file is provided.");
	}

	if (vm.follow("false", "--induce_tree_over_samples") == string("true")) {
		Sankoff sankoff1;
		sankoff1.init(phylo);

 
		if (vm.search("--cost")) {
			sankoff1.load_cost(vm.next(""));
		} else if (vm.search("--cost-identity-aa")) {
			//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
			vector<int> params; // vm["cost-identity-aa"].as<vector<int>>();
			for (int i=0; i<7; i++)
				params.push_back(vm.next(0));
			sankoff1.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
		} else if (vm.search("--cost-identity-dna")) {
			//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
			vector<int> params ; // vm["cost-identity-dna"].as<vector<int>>();
			for (int i=0; i<7; i++)
				params.push_back(vm.next(0));
			sankoff1.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
		} else {
			throw runtime_error("No cost matrix is provided.");
		}
		sankoff1.load_dna(vm.follow("", "--aln"));

		//to use it, map should be rebuild after changing the root
		SubtreeExtractorOverSamples<INode, ValidSampler> subtreeExtractorOverSamples;
		INode phylo_induced = subtreeExtractorOverSamples.run(phylo);
		cerr << "Removed invalid nodes int=" << subtreeExtractorOverSamples.removed_internal_count << " sample=" << subtreeExtractorOverSamples.removed_sample_count << " new tree: size=" << phylo_induced.size << " samples=" << phylo_induced.sample_size << " h=" << phylo_induced.height << endl;

		phylo = phylo_induced;
	}

	//Initialise and run sankoff to assign the best sequence to the internal nodes of the tree
	Sankoff sankoff;
	sankoff.init(phylo);
	
	//TODO merge code with the previous one
	if (vm.search("--cost")) {
		sankoff.load_cost(vm.next(""));
	} else if (vm.search("--cost-identity-aa")) {
		//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
		vector<int> params; // vm["cost-identity-aa"].as<vector<int>>();
		for (int i=0; i<7; i++)
			params.push_back(vm.next(0));
		sankoff.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	} else if (vm.search("--cost-identity-dna")) {
		//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
		vector<int> params ; // vm["cost-identity-dna"].as<vector<int>>();
		for (int i=0; i<7; i++)
			params.push_back(vm.next(0));
		sankoff.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	} else {
		throw runtime_error("No cost matrix is provided.");
	}
	sankoff.load_dna(vm.follow("", "--aln"));

	if (vm.search("--omit-leaf-mutations")) {
		sankoff.omit_leaf_mutations = true;
	}

	sankoff.sankoff();
	cerr << "sankoff done" << endl;

	//if (vm.count("print-alignment")) {
	//	sankoff.print_alignment(vm["print-alignment"].as<string>(), phylo);
	//}
	//string output = vm["out"].as<string>();
	//sankoff.print_mutation_samples(output, phylo);
	//
	//if (vm.count("out-dep") > 0) {
	//	sankoff.print_mutation_dependency(vm["out-dep"].as<string>(), phylo);
	//}

	// //Initialise Mutation dependency calculater using options & run
	// vector<int> important_columns = find_important_columns(vm["min-second-var-count"].as<int>(), sankoff.label_node, sankoff.seq_len, sankoff.char_model);

	if (vm.search("--ilabel")) {
		InternalNodeCustomLabeler<INode> lab(vm.next("inode"));
		lab.label(phylo);
	}

	if (vm.search("--out-as")) {
		sankoff.print_alignment(vm.next(""), phylo);
	}

	if (vm.search("--out-tree")) {
		NodePrinterGeneral<INode> npg;
		ofstream fo(vm.next(""));
		npg.print(fo, phylo);
	}

	// tick_count t1 = tick_count::now();
	cerr << "Done in " << 0 << " seconds" << endl;

	return 0;
}
