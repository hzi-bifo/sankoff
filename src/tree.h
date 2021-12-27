#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <algorithm>
#include <cassert>
#include <regex>
#include <list>
using namespace std;

template<typename T>
ostream& operator<<(ostream& os, const vector<T>& v) {
	os << "[";
	for (auto const &i: v)
		os << i << " ";
	return os << "]";
}

template<typename T, size_t L>
ostream& operator<<(ostream& os, const array<T, L>& v) {
	os << "[";
	for (auto const &i: v)
		os << i << " ";
	return os << "]";
}

template<typename T, typename S>
ostream& operator<<(ostream& os, const map<T,S>& v) {
	os << "[";
	for (auto const &i: v)
		os << i.first <<":" << i.second << " ";
	return os << "]";
}

bool endswith(const string& s, const string& e) {
	if (s.size() < e.size()) return false;
	for (size_t i = s.size() - e.size(), j=0; j < e.size(); i++, j++) {
		if (s[i] != e[j]) return false;
	}
	return true;
}

bool iequals(const string& a, const string& b) {
	for (size_t i = 0, j = 0; i < a.size() || j < b.size(); ) {
		if (i < a.size() && isspace(a[i])) {
			i++;
		} else if (j < b.size() && isspace(b[j])) {
			j++;
		} else if (i < a.size() && j < b.size() && tolower(a[i]) == tolower(b[j])) {
			i++;
			j++;
		} else {
			return false;
		}
	}
	return true;
}

std::string trim(const std::string& str,
                 const std::string& whitespace = " \t\n\r")
{
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    string r = str.substr(strBegin, strRange);
    //cerr << "trim " << r << endl;
    return r;
}

vector<string> split(const string& l, char splitter = '\t') {
	vector<string> x;
	std::istringstream iss(l);
	for (string token; getline(iss, token, splitter); ) {
		x.push_back(token);
	}
	return x;
}


istream& seek_spaces(istream& is, char space = ' ') {
	while (is.peek() == space)
		is.get();
	return is;
}

template<typename DATA=int>
struct Node {
	typedef DATA DataType;
	string label;
	float branch_length;
	string annotation;
	list<Node<DATA>> children;
	int size, height, sample_size;
	int edge_label;//used in jplace

	DATA data;

	Node(string _label="", float _branch_length=0, string _annotation="", 
		const list<Node<DATA>>& _children=list<Node<DATA>>(), 
		int _size=1, int _height=1, int _sample_size=0, int edge_label = -1) :
		label(_label), branch_length(_branch_length), annotation(_annotation), children(_children),
		size(_size), height(_height), sample_size(_sample_size), edge_label(edge_label) {
	}

	void build_tree(istream& is) {
		//cerr << "BT " << (char) is.peek() << endl;
		if (is.peek() == '(') {
			do {
				is.get(); // either first ( or next ,s
				children.push_back(Node<DATA>());
				children.back().build_tree(is);
			} while (is.peek() == ',');
			assert(is.get() == ')');
		}
		label = read_name(is);
		annotation = read_annotation(is);
		branch_length = read_branch_length(is);
		edge_label = read_edge_label(is);
		//
		//if (label == "EPI_ISL_1142103") {
		//	cerr << "Load " << label << " " << branch_length << " " << annotation << endl;
		//}

		sample_size = isLeaf() ? 1 : 0;
		size = height = 1;
		for (auto &i : children) {
			size += i.size;
			sample_size += i.sample_size;
			height = max(height, i.height + 1);
		}
	}

	string read_name(istream& is) {
		string r;
		int c = is.get();
		while (c != ':' && c != ',' && c != '(' && c != ')' && c != '[' && c != ']' && c != ';' && !is.eof()) {
			r += (char) c;
			c = is.get();
		}
		is.putback(c);
		//cerr << " name: " << r << endl;
		return r;
	}

	int read_edge_label(istream& is) {
		if (is.peek() == '{') {
			assert(is.get() == '{');
			int e;
			is >> e;
			assert(is.get() == '}');
			return e;
		}
		return -1;
	}

	float read_branch_length(istream& is) {
		if (is.peek() == ':') {
			assert(is.get() == ':');
			float bl;
			is >> bl;
			//cerr << " bl: " << bl << endl;
			return bl;
		}
		return 0;
	}

	string read_annotation(istream& is) {
		string rr;
		if (is.peek() == '[') {
			char open_brace = is.get();
			if (is.peek() == '&') {
				string r;
				assert(is.get() == '&');
				while (is.peek() != ']')
					r += (char) is.get();
				assert(is.get() == ']');

				//cerr << "annotation: " << label << " : " << r << endl;
				istringstream r_iss(r);
				for (string nv, n, v, location; getline(r_iss, nv, ','); ) {
					istringstream nv_iss(nv);
					getline(nv_iss, n, '=');
					if (n == "location") {
						nv_iss >> location;
						//cerr << "R loc " << label << " " << location << endl;
					} else {
						if (rr.size() > 0) rr += ",";
						rr += nv;
					}
				}
				
			} else {
				is.putback(open_brace);
			}
		}
		return rr;
	}

	vector<string> leafLabels() const {
		vector<string> r;
		if (isLeaf()) {
			r.push_back(label);
			return r;
		}
		for (auto const& n : children) {
			vector<string> rr = n.leafLabels();
			r.insert(r.end(), rr.begin(), rr.end());
		}
		return r;
	}

	bool isLeaf() const {
		return children.size() == 0;
	}

	void name_internal_nodes(int& index) {
		if (!isLeaf()) 
			label = "inode" + to_string(index++);
		for (auto &n : children)
			n.name_internal_nodes(index);
	}

	template<typename T>
	void apply(T& f) {
		f(*this);
		for (auto & c : children) {
			c.apply(f);
		}
	}

	template<typename T>
	void apply_with_info(T& f, Node& parent) {
		f(*this, parent);
		for (auto & c : children) {
			c.apply_with_info(f, *this);
		}
	}

	void update_stat() {
		size = 1;
		height = 1;
		sample_size = isLeaf();
		for (const auto& c : children) {
			size += c.size;
			height = max(height, c.height + 1);
			sample_size += c.sample_size;
		}
	}
};


string get_name(string s) {
	vector<string> x;
	std::istringstream iss(s);
	for (string token; getline(iss, token, '/'); ) {
		x.push_back(token);
	}
	if (x.size() == 4) {
		x.erase(x.begin());
	}
	string r = "";
	for (auto i : x) {
		r += i;
	}
	return r;
}

template<typename NODE>
NODE load_tree(istream& fi) {
	//string l;
	//getline(fi, l);
	//istringstream is(l);
	NODE n;
	n.build_tree(fi);
	cerr << "tree with " << n.size << " nodes and " << n.height << " height loaded " << endl;
	return n;
}

template<typename NODE>
struct node_rename_by_map {
	const map<string, string>& label_map;

	node_rename_by_map(const map<string, string>& _label_map) : label_map(_label_map) {}

	void operator()(NODE& n) {
		if (n.isLeaf()) {
			if(label_map.find(n.label) == label_map.end()) {
				cerr << "Label " << n.label << " not present in map" << endl;
			}
			assert(label_map.find(n.label) != label_map.end());
			//cerr << " L " << n.label << " -> " << label_map.find(n.label)->second << endl;
			n.label = label_map.find(n.label)->second;
		}
	}
};

template<typename NODE>
NODE load_nexus_tree(ifstream& fi) {
	//cerr << "  nexus: begin " << endl;
	map<string, string> id_to_name;
	for (string s; getline(fi, s); ) {
		//cerr << "  nexus: read " << s.substr(0, 40) << endl;
		if (iequals(s, "begin trees;")) {
			fi >> s;
			//cerr << " after beg tree " << s << endl;
			if (iequals(trim(s), "translate")) {
				while (trim(s) != ";") {
					vector<string> row = split(trim(s), '\t');
					if (row.size() != 2)
						row = split(trim(s), ' ');
					if (row.size() == 2) {
						string name = trim(row[1]);
						if (name.size() > 0 && name[name.size()-1] == ',')
							name = name.substr(0, name.size()-1);
						id_to_name[row[0]] = name;
					}
					getline(fi, s);
				}
				fi >> s;
			}
			if (iequals(s, "tree")) {
				for (int c, i=0; (c=fi.get()) != '='; i++) {
					//cerr << "C " << (char) c << endl;
					assert(i < 50);
				}
				seek_spaces(fi);
				if (fi.peek() == '[') {
					while (fi.get() != ']')
						;
				}
				seek_spaces(fi);
				assert(fi.peek() == '(');
				break;
			} else {
				assert(1 != 1);
			}
		} else {
			//before begin trees
			//cerr << "LN bt != " << s.substr(0, 15) << endl;
		}
	}
	NODE n = load_tree<NODE>(fi);
	if (id_to_name.size() > 0) {
		node_rename_by_map<NODE> renamer(id_to_name);
		n.apply(renamer);
	}
	return n;
}

template<typename NODE>
NODE load_tree(string fn) {
	//cerr << "loading " << fn << " ..." << endl;
	ifstream fi(fn);
	NODE n = endswith(fn, ".nexus") ? load_nexus_tree<NODE>(fi) : load_tree<NODE>(fi);
	return n;
}

int indexOf(const vector<string> v, const std::initializer_list<string>& keys) {
	for (auto const &k : keys) {
		if ( find(v.begin(), v.end(), k) != v.end())
			return find(v.begin(), v.end(), k) - v.begin();
	}
	cerr << "Warning: not found indexOf "  << v << " " << vector<string>(keys.begin(), keys.end()) << endl;
	//assert(1 != 1);
	return -1;
}

bool startsWith(string s, string w) {
	//cerr << "SW " << s << " " << w << endl;
	return s.find(w) == 0;
}

template<typename NODE>
struct NodePrinterAbstractClass {
	bool force_print_location;
	bool print_internal_node_labels, allow_print_annotation;

	NodePrinterAbstractClass(bool force_print_location, bool print_internal_node_labels, bool allow_print_annotation = true) : 
		force_print_location(force_print_location), print_internal_node_labels(print_internal_node_labels), allow_print_annotation(allow_print_annotation) {
	}

	ostream& print_node_info(ostream&os, const NODE& n) const {
		if (n.label != "" && (print_internal_node_labels || (n.isLeaf()))) {
			os << n.label;
		}
		if (allow_print_annotation && (n.annotation != "" || (n.location != NODE::StateType::unknown && force_print_location) /*|| sankoff_value[0] > -1*/)) {
			os << "[&";
			bool empty = true;
			if (n.annotation != "") {
				os << n.annotation;
				empty = false;
			}
			if (n.location != NODE::StateType::unknown && force_print_location) {
				if (!empty)
					os << ",";
				os << "location=" << n.location;
				empty = false;
			}
			//os << "," << "s=" << sankoff_value[0] << "|" << sankoff_value[1];
			os << "]";
		}
		return os << ":" << n.branch_length;
	}

};

template<typename NODE>
struct NodePrinterGeneral : public NodePrinterAbstractClass<NODE> {

	NodePrinterGeneral(bool _force_print_location = true, bool _print_internal_node_labels = true, bool allow_print_annotation = true) : 
		NodePrinterAbstractClass<NODE>(_force_print_location, _print_internal_node_labels, allow_print_annotation) {}

	ostream& print(ostream& os, const NODE& n) const {
		if (n.children.size() > 0) {
			os << "(";
			bool first = true;
			for (auto &c : n.children) {
				if (!first)
					os << ",";
				first = false;
				print(os, c);
			}
			os << ")";
		}
		this->print_node_info(os, n);
		return os;
	}

};

template<typename NODE>
struct SamplePrinter {
	int printed_count = 0;
	SamplePrinter() {}

	void run(const NODE& n, string fn) {
		//cerr << "SamplePrinter::run() " << fn << " " << n.label << endl;
		ofstream fo(fn);
		//cerr << "SamplePrinter::run() file created " << endl;
		dfs(n, fo);
	}

	void dfs(const NODE& n, ofstream& fo) {
		//cerr << "sample printer " << n.label << endl;
		if (n.isLeaf()) {
			printed_count++;
			fo << n.label << endl;
		}
		for (auto const& c : n.children) {
			dfs(c, fo);
		}
	}
};

template<typename NODE>
struct InternalNodeLabeler {
	int cnt;
	InternalNodeLabeler() : cnt(0) {}
	void run(NODE& n) {
		cnt=0;
		label(n);
	}

	void label(NODE&n) {
		if (!n.isLeaf()) {
			n.label = "inode" + to_string(cnt++);
		}
		for (auto &c : n.children) {
			label(c);
		}
	}
};

template<typename NODE, typename SampleIncluder>
struct SubtreeExtractorOverSamples {

	int removed_internal_count, removed_sample_count;
	SampleIncluder sampleIncluder;

	SubtreeExtractorOverSamples(const SampleIncluder& sampleIncluder = SampleIncluder()) : removed_internal_count(0), removed_sample_count(0), 
		sampleIncluder(sampleIncluder) {
	}
	
	NODE run(const NODE& n) {
		NODE r = sample(n).first;
		if (r.children.size() == 1) 
			return r.children.back();
		return r;
	}

	//rebuild subtree basd on sampled/non-sampled, each internal node should have at least two children.
	pair<NODE, bool> sample(const NODE& n) {
		//height, size should be recalculated
		NODE r(n.label, n.branch_length, n.annotation, list<NODE>(), 1, 1, n.isLeaf() ? 1 : 0);
		for (auto &c: n.children) {
			pair<NODE, bool> c_c_ = sample(c);
			if (c_c_.second == false) {
				if (c.isLeaf())
					removed_sample_count++;
				else
					removed_internal_count++;
				continue;
			}
			NODE& c_c = c_c_.first;
			if (c_c.children.size() == 1 && !c.isLeaf()) {
				//cerr << "node removed " << c.label << " " << c.branch_length << endl;
				removed_internal_count++;
				for (auto& cc: c_c.children) {
					cc.branch_length += c.branch_length;
					r.height = max(r.height, cc.height + 1);
					r.size += cc.size;
					//r.children.push_back(move(cc));
					r.sample_size += cc.sample_size;
					r.children.push_back(cc);
				}
			} else {
				//normal:
				r.height = max(r.height, c_c.height + 1);
				r.size += c_c.size;
				r.sample_size += c_c.sample_size;
				//r.children.push_back(move(c_c));
				r.children.push_back(c_c);
			}
		}
		//return make_pair(r, r.children.size() > 0 || sample_names.find(n.label) != sample_names.end());
		return make_pair(r, r.children.size() > 0 || sampleIncluder(n));
	}
};

template<typename NODE>
struct SingleChildInternalNodeRemover {

	int removed_internal_count;

	SingleChildInternalNodeRemover() : removed_internal_count(0) {
	}
	
	NODE run(const NODE& n) {
		removed_internal_count = 0;
		return dfs(n);
	}

	NODE dfs(const NODE& n) {
		NODE r(n.label, n.branch_length, n.annotation, list<NODE>(), n.location, 1, 1, n.isLeaf() ? 1 : 0);
		for (auto &c: n.children) {
			NODE c_c = dfs(c);
			//normal:
			r.height = max(r.height, c_c.height + 1);
			r.size += c_c.size;
			r.sample_size += c_c.sample_size;
			//r.children.push_back(move(c_c));
			r.children.push_back(c_c);
		}
		if (r.children.size() == 1) {
			removed_internal_count++;
			r.children.begin()->branch_length += r.branch_length;
			//r = std::move(*r.children.begin());
			//TODO: we should make it faster!
			//cerr << r.children.begin()->label << endl;
			NODE tmp = *r.children.begin();
			r = tmp;
		}
		return r;
	}
};

template<typename NODE, typename DFS_ACTION>
struct TreeDFSGeneral {

	DFS_ACTION action;
	TreeDFSGeneral(DFS_ACTION action) : action(action) {}

	void dfs(NODE& n) {
		action.visit(n);

		for (auto & c : n.children) {
			dfs(c);
		}
		action.finish(n);
	}
};


template<typename NODE>
struct NodePrinterNexus {

	vector<pair<int, string>> index_name;
	int current_index;

	NodePrinterNexus() {}

	void print_node_info(ostream& os, const NODE& n) {
		if (n.label == "EPI_ISL_753872") 
			cerr << "going to pring EPI_ISL_753872 ++ " << endl;

		if (n.label != "") { //  && n.isLeaf()
			int l = current_index++;
			index_name.push_back(make_pair(l, n.label));
			os << l;
		}
		if (n.annotation != "") {
			os << "[&";
			bool empty = true;
			if (n.annotation != "") {
				os << n.annotation;
				empty = false;
			}
			if (n.location != NODE::StateType::unknown) {
				if (!empty)
					os << ",";
				os << "location=" << n.location;
				empty = false;
			}
			//os << "," << "s=" << sankoff_value[0] << "|" << sankoff_value[1];
			os << "]";
		}
		os << ":" << n.branch_length;
	}

	ostream& print(ostream& os, const NODE& n) {
			if (n.label == "EPI_ISL_753872") 
				cerr << "going to pring EPI_ISL_753872 " << endl;

		if (n.children.size() > 0) {
			os << "(";
			bool first = true;
			for (auto &c : n.children) {
				if (!first)
					os << ",";
				first = false;
				print(os, c);
			}
			os << ")";
		}
			if (n.label == "EPI_ISL_753872") 
				cerr << "going to pring EPI_ISL_753872 + " << endl;
		print_node_info(os, n);
		return os;
	}

	void run(ostream& os, const NODE& n) {
		index_name.clear();
		current_index = 1;
		print(os, n);
	}

};


template<typename NODE>
void save_nexus_tree(ostream& os, const NODE& n) {
	//cerr << "  nexus: begin " << endl;
	ostringstream tree_oss;
	NodePrinterNexus<NODE> npn;
	npn.run(tree_oss, n);

	os << "#NEXUS" << endl << endl
	   << "Begin taxa;" << endl
	   << "\tDimensions ntax=" << npn.index_name.size() << ";" << endl
	   << "\tTaxlabels" << endl;
	for (auto & iname : npn.index_name) {
		os << "\t\t" << iname.second << endl;
	}
	os << "\t\t;" << endl << "End;" << endl << endl
	   << "Begin trees;" << endl
	   << "\tTranslate";
	bool first = true;
	for (auto & iname : npn.index_name) {
		if (!first)
			os << ",";
		os << endl
		   << "\t\t" << iname.first << " " << iname.second;
		first = false;
	}
	os << endl << "\t;" << endl;

	os << "tree TREE1 = [&R] " 
	   << tree_oss.str() << ";" << endl
	   << "End;" << endl;
}

