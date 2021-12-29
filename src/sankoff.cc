#include "tree.h"
#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb.h>
#include <boost/program_options.hpp>

using namespace oneapi::tbb;
namespace po = boost::program_options;

#include "sankoff.h"


int main(int argc, char* argv[]) {
	
	po::options_description desc("Fast Sankoff algorithm with parallel options.");
	desc.add_options()
	    ("help", "print help")
	    ("tree", po::value<string>(), "input tree file.")
	    ("nexus", po::value<string>(), "input nexus file.")
	    ("aln", po::value<string>()->required(), "aligned sequence file")
	    ("cost", po::value<string>(), "cost function file name")
	    ("cost-identity-aa", po::value<vector<int>>(), "Cost matrix is for chars of amino acids in addition to X and -. The parameters are cost of A-B where A and B are either amino acids (aa) X or gap with following order [identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap].")
	    ("cost-identity-dna", po::value<vector<int>>(), "Cost matrix is for chars of nucleic acid in addition to X and -. The parameters are cost of A-B where A and B are either nucleid acids (na) X or gap with following order [identical na] [non-identical na] [na-X] [X-X] [gap-X] [na-gap] [gap-gap].")
		//TODO: some default and simple cost functions
//	    ("print-alignment", po::value<string>(), "print alignment")
	    ("omit-leaf-mutations", "omit mutations happen at leaf nodes")
	    ("nthread", po::value<int>()->default_value(oneapi::tbb::this_task_arena::max_concurrency()), "change number of default threads")
	    ("induce_tree_over_samples", po::value<bool>()->default_value(true), "remove nodes not in the alignment")
	    ("asf", po::value<string>(), "print sequences of ancestral and leaf nodes to this file.")
	    ("ilabel", po::value<string>()->default_value("inode"), "Assign label to internal nodes. The argument is the prefix.")
	;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);

	if (vm.count("help"))
	{
		cout << "Usage: " << argv[0] << " [options] ...\n";
		cout << desc;
		return 0;
	}
	po::notify(vm);    

	oneapi::tbb::global_control global_limit(oneapi::tbb::global_control::max_allowed_parallelism, vm["nthread"].as<int>());


	tick_count t0 = tick_count::now();

	//VM vm(argv);

	//Build tree from MSA removing excess nodes
	INode phylo;
	if (vm.count("tree") > 0) {
		ifstream fi(vm["tree"].as<string>());
		phylo = load_tree<INode>(fi);
	} else 	if (vm.count("nexus") > 0) {
		ifstream fi(vm["nexus"].as<string>());
		phylo = load_nexus_tree<INode>(fi);
	} else {
		throw runtime_error("No tree file is provided.");
	}

	if (vm["induce_tree_over_samples"].as<bool>()) {
		Sankoff sankoff1;
		sankoff1.init(phylo);

 
		if (vm.count("cost")) {
			sankoff1.load_cost(vm["cost"].as<string>());
		} else if (vm.count("cost-identity-aa")) {
			//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
			vector<int> params = vm["cost-identity-aa"].as<vector<int>>();
			sankoff1.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
		} else if (vm.count("cost-identity-dna")) {
			//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
			vector<int> params = vm["cost-identity-dna"].as<vector<int>>();
			sankoff1.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
		} else {
			throw runtime_error("No cost matrix is provided.");
		}
		sankoff1.load_dna(vm["aln"].as<string>());

		//to use it, map should be rebuild after changing the root
		SubtreeExtractorOverSamples<INode, ValidSampler> subtreeExtractorOverSamples;
		INode phylo_induced = subtreeExtractorOverSamples.run(phylo);
		cerr << "Removed invalid nodes int=" << subtreeExtractorOverSamples.removed_internal_count << " sample=" << subtreeExtractorOverSamples.removed_sample_count << " new tree: size=" << phylo_induced.size << " samples=" << phylo_induced.sample_size << " h=" << phylo_induced.height << endl;

		phylo = phylo_induced;
	}

	//Initialise and run sankoff to assign the best sequence to the internal nodes of the tree
	Sankoff sankoff;
	sankoff.init(phylo);
	

	if (vm.count("cost")) {
		sankoff.load_cost(vm["cost"].as<string>());
	} else if (vm.count("cost-identity-aa")) {
		//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
		vector<int> params = vm["cost-identity-aa"].as<vector<int>>();
		sankoff.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	} else if (vm.count("cost-identity-dna")) {
		//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
		vector<int> params = vm["cost-identity-dna"].as<vector<int>>();
		sankoff.char_model.init_aa(params[0], params[1], params[2], params[3], params[4], params[5], params[6]);
	} else {
		throw runtime_error("No cost matrix is provided.");
	}
	sankoff.load_dna(vm["aln"].as<string>());

	if (vm.count("omit-leaf-mutations")) {
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


	if (vm.count("ilabel")) {
		InternalNodeCustomLabeler<INode> lab(vm["ilabel"].as<string>());
		lab.label(phylo);
	}

	if (vm.count("out-as")) {
		sankoff.print_alignment(vm["out-as"].as<string>(), phylo);
	}

	if (vm.count("out-tree")) {
		NodePrinterGeneral<INode> npg;
		ofstream fo(vm["out-tree"].as<string>());
		npg.print(fo, phylo);
	}

	tick_count t1 = tick_count::now();
	cerr << "Done in " << (t1-t0).seconds() << " seconds" << endl;

	return 0;
}
