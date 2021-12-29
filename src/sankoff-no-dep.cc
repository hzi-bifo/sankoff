#include "tree.h"

struct CharModel {
	int na_count;

	int na_to_int(char na) const {
		na = toupper(na);
		if (!__char_is_valid[(int)na]) {
			cerr << "Invalid char '" << na << "'" << endl;
		}
		assert(__char_is_valid[(int)na] == true);
		return __char_to_int[(int)na];
	}

	int na_to_int(string na) {
		return na_to_int(na[0]);
	}

	char int_to_char(int na) {
		assert(na >= 0 && na < na_count);
		return __int_to_char[na];
	}

	static const int MAX_NA_COUNT = 30;

	int cost[MAX_NA_COUNT][MAX_NA_COUNT];
	int __char_to_int[255];
	char __int_to_char[MAX_NA_COUNT];
	bool __char_is_valid[255];

	//TODO check cost file error including:
	//  1) names should be 1 char 
	//  2) row names and column names should be equal (no new row and equal number)
	void load_cost(string fn) {
		ifstream f(fn);
		vector<string> names;
		for (size_t i=0; i<255; i++)
			__char_is_valid[i] = false;
		for (string line, v; getline(f, line); ) {
			istringstream is(line.c_str());
			if (names.size() == 0) {
				for (string s; is >> s; )
					names.push_back(s);
				na_count = names.size();
				assert(na_count < MAX_NA_COUNT);
				for (size_t i=0; i<names.size(); i++) {
					__int_to_char[i] = toupper(names[i][0]);
					__char_to_int[toupper(names[i][0])] = i;
					__char_is_valid[(int)toupper(names[i][0])] = true;
				}
			} else {
				is >> v;
				for (int vv, i=0; is >> vv; i++) {
					cost[na_to_int(v)][na_to_int(names[i])] = vv;
				}
			}
		}
	}

	void init_by_vector(vector<string> main_names, string X, string gap, int identical_cost, int non_identical_cost, int X_aa, int XX, int Xgap, int gap_aa, int gapgap) {
		vector<string> names = main_names;
		names.push_back(X);
		names.push_back(gap);
		na_count = (int) names.size();
		for (size_t i=0; i<255; i++)
			__char_is_valid[i] = false;
		for (size_t i = 0 ; i < names.size(); i++) {
			__int_to_char[i] = toupper(names[i][0]);
			__char_to_int[toupper(names[i][0])] = i;
			__char_is_valid[(int)toupper(names[i][0])] = true;
		}
		for (string a : names) {
			for (string b : names) {
				int var = -1;
				if (a == X || a == gap || b == X || b == gap) {
					if (a == X) {
						if (b == X) {
							var = XX;
						} else if (b == gap) {
							var = Xgap;
						} else {
							var = X_aa;
						}
					} else if (a == gap) {
						if (b == X) {
							var = Xgap;
						} else if (b == gap) {
							var = gapgap;
						} else {
							var = gap_aa;
						}
					}
				} else {
					if (a == b) {
						var = identical_cost;
					} else {
						var = non_identical_cost;
					}
				}
				cost[na_to_int(a)][na_to_int(b)] = var;
			}

		}
	}
	//[identical aa] [non-identical aa] [aa-X] [X-X] [gap-X] [aa-gap] [gap-gap]
	void init_aa(int identical_cost, int non_identical_cost, int X_aa, int XX, int Xgap, int gap_aa, int gapgap) {
		vector<string> names{"A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"};
		init_by_vector(names, "X", "-", identical_cost, non_identical_cost, X_aa, XX, Xgap, gap_aa, gapgap);
	}

	void init_na(int identical_cost, int non_identical_cost, int X_aa, int XX, int Xgap, int gap_aa, int gapgap) {
		vector<string> names{"A", "T", "C", "G"};
		init_by_vector(names, "X", "-", identical_cost, non_identical_cost, X_aa, XX, Xgap, gap_aa, gapgap);
	}

};

struct Data {
	bool valid;
	//int sankoff_value[CharModel::MAX_NA_COUNT];
	vector<char> seq;
	size_t node_index;
	void load(string s) {
		seq = vector<char>(s.begin(), s.end());
	}
};
typedef Node<Data> INode;

struct SankoffParallelTaskForLocation {

	INode& root;
	//sankoff_value[node_index][na]
	vector<array<int,CharModel::MAX_NA_COUNT>> sankoff_value;
	CharModel& char_model;


	SankoffParallelTaskForLocation(INode& root, CharModel& char_model) : root(root), char_model(char_model) {
		sankoff_value.resize(root.size);
	}

	void run(int loc) {
		sankoff1(root, loc);
		sankoff2(root, -1, loc);
	}

	void sankoff1(INode& n, int loc) {
		//cerr << "s1 " << n.label << " " << n.data.valid << " " << n.data.seq << " " << endl;
		if (n.children.size() == 0) {
			for (int l = 0; l < char_model.na_count; l++) {
				//n.data.sankoff_value[l] = 1<<25;
				sankoff_value[n.data.node_index][l] = 1<<25;
			}
			sankoff_value[n.data.node_index][char_model.na_to_int(n.data.seq[loc])] = 0;
		} else {
			for (auto&c : n.children) {
				if (c.data.valid) {
					sankoff1(c, loc);
				}
			}
			for (int l = 0; l<char_model.na_count; l++) {
				int v = 0;
				for (auto& c : n.children) {
					if (c.data.valid) {
						int v_with_child = 1<<25;
						for (int l2=0; l2<char_model.na_count; l2++)
							//v_with_child = min(v_with_child, char_model.cost[l][l2] + c.data.sankoff_value[l2]);
							v_with_child = min(v_with_child, char_model.cost[l][l2] + sankoff_value[c.data.node_index][l2]);
						v += v_with_child;
					}
				}
				sankoff_value[n.data.node_index][l] = v;
			}
		}
		/*
		   cerr << "sankoff1 " << n.label << " ";
		for (int i=0; i<NA_COUNT; i++)
			cerr << n.data.sankoff_value[i] << " ";
		cerr << endl;
		*/
	}

	void sankoff2(INode& n, int l, int loc) {
		if (l == -1) {
			int min_l = 0; //default is NON_GERMANY
			for (int l = 0; l<char_model.na_count; l++)
				//if (n.data.sankoff_value[l] < n.data.sankoff_value[min_l])
				if (sankoff_value[n.data.node_index][l] < sankoff_value[n.data.node_index][min_l])
					min_l = l;
			sankoff2(n, min_l, loc);
		} else {
			//if ((int)n.data.seq.size() != seq_len) {
			//	n.data.seq.resize(seq_len);
			//}
			n.data.seq[loc] = char_model.int_to_char(l);
			//cerr << "A " << n.label << " loc=" << loc << " =" << l << endl;
			for (auto&c : n.children) {
				if (c.data.valid) {
					//double v_with_child = 1e10;
					int min_l2 = l;
					for (int l2=0; l2<char_model.na_count; l2++) {
						//float cost_new = char_model.cost[l][l2] + c.data.sankoff_value[l2], 
						float cost_new = char_model.cost[l][l2] + sankoff_value[c.data.node_index][l2], 
							cost_old = char_model.cost[l][min_l2] + sankoff_value[c.data.node_index][min_l2];
						//if (cost[l][l2] + n.sankoff_value[l2] < cost[l][min_l2] + n.sankoff_value[min_l2])
						if (cost_new < cost_old || (cost_new == cost_old && l2 == l))
							min_l2 = l2;
						//v_with_child = min(v_with_child, cost[l][l2] + n.sankoff_value[l2]);
					}
					sankoff2(c, min_l2, loc);
				}
			}
		}
	}


};

struct SankoffParallelTaskForLocationExecutor {

	INode& root;
	//sankoff_value[node_index][na]
	vector<array<int,CharModel::MAX_NA_COUNT>> sankoff_value;
	CharModel& char_model;

	SankoffParallelTaskForLocationExecutor(INode& root, CharModel& char_model) : root(root), char_model(char_model) {
		sankoff_value.resize(root.size);
	}

	void operator()( const blocked_range<size_t>& r ) const {
		SankoffParallelTaskForLocation sankoff_task(root, char_model);
		for( size_t i=r.begin(); i!=r.end(); ++i ) {
			sankoff_task.run(i);
			cerr << "\rS loc=" << i << "/" << (r.end() - r.begin());
		}
	}

};


struct Sankoff {
	map<string, INode*> label_node;
	INode* root;
	int seq_len;
	CharModel char_model;
	bool omit_leaf_mutations = false;

	void load_cost(string fn) {
		char_model.load_cost(fn);
	}

	void init(INode& _root) {
		root = &_root;
		int node_index = 0;
		build_map(*root, node_index);
		cerr << "build_map done with " << label_node.size() << " samples" << " node_index=" << node_index << endl;
		assert(node_index == _root.size);

	}

	void build_map(INode& n, int &node_index) {
		//cerr << "build_map " << n.label << " " << label_node.size() << endl;
		n.data.node_index = node_index++;
		n.data.valid = true;
		if (n.isLeaf()) {
			label_node[n.label] = &n;
		}
		for (auto & c : n.children) {
			build_map(c, node_index);
		}
	}

	void load_dna(string fn) {
		cerr << "loading sequences for " << label_node.size() << " nodes" << endl;
		assert(label_node.size() != 0);
		//for (auto & m : label_node)
		//	cerr << m.first << " ";
		//cerr << endl;

		//seq_len = 40000;
		//for (auto & i : label_node) {
		//	i.second->data.seq = vector<char>(seq_len, 'A');
		//}
		//return;

		seq_len = -1;
		/*
		ifstream file(fn, ios_base::in | ios_base::binary);
		cerr << "file opened " << fn << endl;
		io::filtering_streambuf<io::input> in;
		//in.push(gzip_decompressor());
		in.push(io::lzma_decompressor());
		in.push(file);
		std::istream incoming(&in);
		istream& file_alignment = endswith(fn,".xz") ? incoming : file; 
		*/
		ifstream file_alignment(fn);

		string seq, id;
		for (string line; getline(file_alignment, line); ) {
			//cerr << "L " << line << endl;
			if (line.size() == 0) continue;
			if (line[0] == '>') {
				if (id != "" && label_node.find(id) != label_node.end()) {
					//cerr << "S " << id << endl;
					//cout << ">" << id << endl;
					//cout << seq << endl;
					//samples_found.insert(id);
					label_node[id]->data.load(seq);
					//cerr << "load seq " << seq << " for node " << id << endl;
					if (seq_len == -1) seq_len = seq.size();
					assert(seq_len == (int) seq.size());
				} else {
					cerr << "ID not in tree " << id << endl;
				}
				//id = split(line.substr(1), '|')[0];
				id = line.substr(1);
				seq = "";
				//cerr << "new id (" << id << ") " << id.size() << endl;
			} else {
				seq += line;
			}
		}
		if (id != "" && label_node.find(id) != label_node.end()) {
			//cout << ">" << id << endl;
			//cout << seq << endl;
			//samples_found.insert(id);
			label_node[id]->data.load(seq);
		}
		cerr << "sequence files loaded" << endl;

		string leafs_without_seq = "";
		int leafs_without_seq_count = 0, leafs_with_seq_count = 0;
		for (auto & i : label_node) {
			if ((int)i.second->data.seq.size() != seq_len) {
				leafs_without_seq_count++;
				//keeping label of first 10 leaf nodes without sequence for printing error
				if (leafs_without_seq_count < 10) {
					if (leafs_without_seq.size() > 0)
						leafs_without_seq += ", ";
					leafs_without_seq += i.second->label;
				}
				//cerr << "E " << i.second->label << " no sequence data" << endl;
				i.second->data.valid = false;
			} else {
				leafs_with_seq_count++;
			}
			//assert((int)i.second->data.seq.size() == seq_len);
			//for (auto & c : i.second->data.seq) {
			//	c = toupper(c);
			//	//if (c != 'A' && c != 'T' && c != 'C' && c != 'G' && c != '-' && c != 'N') {
			//	//	cerr << "invalid char " << i.second->label << " " << c << endl;
			//	//	assert( 1 != 1 );
			//	//}
			//}
		}
		if (leafs_without_seq_count > 0) {
			cerr << "E no sequence data for " << leafs_without_seq_count << " leaf nodes (i.e. " << leafs_without_seq << ")" << " seq_len = " << seq_len << endl;
			cerr << "  with seq " << leafs_with_seq_count << endl;
		}
	}

	void print_alignment(string fn, const INode& n) {
		ofstream fo(fn);
		cerr << "print_alignment " << fn << endl;
		print_alignment(fo, n);
	}

	void print_alignment(ofstream& fo, const INode& n) {
		string s(n.data.seq.begin(), n.data.seq.end());
		fo << ">" << n.label << endl << s << endl;
		for (auto const& c : n.children)
			print_alignment(fo, c);
	}

//	void sankoff1(INode& n, int loc) {
//		//cerr << "s1 " << n.label << " " << n.data.valid << " " << n.data.seq << " " << endl;
//		if (n.children.size() == 0) {
//			for (int l = 0; l < char_model.NA_COUNT; l++)
//				n.data.sankoff_value[l] = 1<<25;
//			n.data.sankoff_value[char_model.na_to_int(n.data.seq[loc])] = 0;
//		} else {
//			for (auto&c : n.children) {
//				if (c.data.valid) {
//					sankoff1(c, loc);
//				}
//			}
//			for (int l = 0; l<char_model.NA_COUNT; l++) {
//				int v = 0;
//				for (auto& c : n.children) {
//					if (c.data.valid) {
//						int v_with_child = 1<<25;
//						for (int l2=0; l2<char_model.NA_COUNT; l2++)
//							v_with_child = min(v_with_child, char_model.cost[l][l2] + c.data.sankoff_value[l2]);
//						v += v_with_child;
//					}
//				}
//				n.data.sankoff_value[l] = v;
//			}
//		}
//		/*
//		   cerr << "sankoff1 " << n.label << " ";
//		for (int i=0; i<NA_COUNT; i++)
//			cerr << n.data.sankoff_value[i] << " ";
//		cerr << endl;
//		*/
//	}
//
//	void sankoff2(INode& n, int l, int loc) {
//		if (l == -1) {
//			int min_l = 0; //default is NON_GERMANY
//			for (int l = 0; l<char_model.NA_COUNT; l++)
//				if (n.data.sankoff_value[l] < n.data.sankoff_value[min_l])
//					min_l = l;
//			sankoff2(n, min_l, loc);
//		} else {
//			//if ((int)n.data.seq.size() != seq_len) {
//			//	n.data.seq.resize(seq_len);
//			//}
//			n.data.seq[loc] = char_model.int_to_char(l);
//			//cerr << "A " << n.label << " loc=" << loc << " =" << l << endl;
//			for (auto&c : n.children) {
//				if (c.data.valid) {
//					//double v_with_child = 1e10;
//					int min_l2 = l;
//					for (int l2=0; l2<char_model.NA_COUNT; l2++) {
//						float cost_new = char_model.cost[l][l2] + c.data.sankoff_value[l2], 
//							cost_old = char_model.cost[l][min_l2] + c.data.sankoff_value[min_l2];
//						//if (cost[l][l2] + n.sankoff_value[l2] < cost[l][min_l2] + n.sankoff_value[min_l2])
//						if (cost_new < cost_old || (cost_new == cost_old && l2 == l))
//							min_l2 = l2;
//						//v_with_child = min(v_with_child, cost[l][l2] + n.sankoff_value[l2]);
//					}
//					sankoff2(c, min_l2, loc);
//				}
//			}
//		}
//	}

	void assign_memory(INode& n) {
		if ((int)n.data.seq.size() != seq_len)
			n.data.seq.resize(seq_len, 0);
		for (auto& c: n.children)
			assign_memory(c);
	}

	void sankoff_pre(INode& n, int& valid_nodes, int& invalid_nodes) {
		//validity of leaf nodes should be filled when dna is loaded
		if (!n.isLeaf())
			n.data.valid = false;
		for (auto& c: n.children) {
			sankoff_pre(c, valid_nodes, invalid_nodes);
			n.data.valid |= c.data.valid;
		}
		if (n.data.valid) {
			valid_nodes++;
		} else {
			invalid_nodes++;
		}
	}

	void sankoff() {
		sankoff(*root);
	}

	void sankoff(INode& n) {
		//cerr << "assigning memory " << endl;
		assign_memory(n);
		cerr << "sankoff " << n.label << " " << seq_len << endl;
		int valid_nodes = 0, invalid_nodes = 0;
		sankoff_pre(n, valid_nodes, invalid_nodes);
		cerr << "Tree valid=" <<valid_nodes << " " << "invalid=" << invalid_nodes << endl;
		cerr << "S ...";
		/*
		for (int i=0; i<seq_len; i++) {
			sankoff1(n, i);
			sankoff2(n, -1, i);
			cerr << "\rS loc=" << i << "/" << seq_len;
			//for (int j=0; j<NA_COUNT; j++)
			//	cerr << n.data.sankoff_value[j] << " ";
			//cerr << endl;

		}
		*/
		parallel_for(blocked_range<size_t>(0, seq_len), SankoffParallelTaskForLocationExecutor(n, char_model));
		cerr << endl;
	}
};


struct ValidSampler {
	bool operator()(const INode& n) const {
		return n.isLeaf() ? n.data.valid : false;
		//return n.data.valid;
	}
};

struct VMValue {
	void* value;
	template<typename T>
	T as() {
		return *static_cast<T*>(value);
	}
	VMValue(void* value=0) : value(value) {}
};

vector<int> find_important_columns(int min_second_variant_count, const map<string, INode*>& label_node, int seq_len, const CharModel& char_model) {
	vector<int> ret;
	for (int i=0; i<seq_len; i++) { //iterate over each position
		int na_cnt[CharModel::MAX_NA_COUNT] = {}; // vector to hold count values for each AA
		for (auto& ln : label_node) {
			na_cnt[char_model.na_to_int(ln.second->data.seq[i])]++;
			//cerr << ln.second->data.seq[i]; //char_model.na_to_int(ln.second->data.seq[i])
		}
		sort(na_cnt, na_cnt+char_model.na_count, greater<int>());// sort counts
		if (na_cnt[1] >= min_second_variant_count) { // check second highest value against threshold value
			//cerr << "Important column " << i << " " << na_cnt[0] << " " << na_cnt[1] << endl;
			ret.push_back(i); // add to important columns
		}
	}
	return ret;
}



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
	
	sankoff.load_cost(vm["cost"].as<string>());
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

	tick_count t1 = tick_count::now();
	cerr << "Done in " << (t1-t0).seconds() << " seconds" << endl;

	return 0;
}
