#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <boost/xpressive/xpressive.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <algorithm>
#include <regex>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>

using namespace std;
using namespace boost::xpressive;

class DifferenceSite {
private:
    string id;
    int in_alignment_site;
    int in_seq_site;
    char char_in_different;

public:
    const string &getId() const {
        return id;
    }

    void setId(const string &id) {
        DifferenceSite::id = id;
    }

    int getIn_alignment_site() const {
        return in_alignment_site;
    }

    void setIn_alignment_site(int in_alignment_site) {
        DifferenceSite::in_alignment_site = in_alignment_site;
    }

    int getIn_seq_site() const {
        return in_seq_site;
    }

    void setIn_seq_site(int in_seq_site) {
        DifferenceSite::in_seq_site = in_seq_site;
    }

    char getChar_in_different() const {
        return char_in_different;
    }

    void setChar_in_different(char char_in_different) {
        DifferenceSite::char_in_different = char_in_different;
    }
};

class Primer {
private:
    string id;
    int start;
    int length;
    string tm;
    string gc;
    string any;
    string three;
    string seq;
    int diff_site_number;
    string diff_3_end_flag;
    vector<DifferenceSite> diff_vector;
    double score;
    int start_in_alignment;
    int end_in_alignment;
    int deletion_site_number;

public:
    const string &getId() const {
        return id;
    }

    void setId(const string &id) {
        Primer::id = id;
    }

    int getStart() const {
        return start;
    }

    void setStart(int start) {
        Primer::start = start;
    }

    int getLength() const {
        return length;
    }

    void setLength(int length) {
        Primer::length = length;
    }

    const string &getTm() const {
        return tm;
    }

    void setTm(const string &tm) {
        Primer::tm = tm;
    }

    const string &getGc() const {
        return gc;
    }

    void setGc(const string &gc) {
        Primer::gc = gc;
    }

    const string &getAny() const {
        return any;
    }

    void setAny(const string &any) {
        Primer::any = any;
    }

    const string &getThree() const {
        return three;
    }

    void setThree(const string &three) {
        Primer::three = three;
    }

    const string &getSeq() const {
        return seq;
    }

    void setSeq(const string &seq) {
        Primer::seq = seq;
    }

    int getDiff_site_number() const {
        return diff_site_number;
    }

    void setDiff_site_number(int diff_site_number) {
        Primer::diff_site_number = diff_site_number;
    }

    const string &getDiff_3_end_flag() const {
        return diff_3_end_flag;
    }

    void setDiff_3_end_flag(const string &diff_3_end_flag) {
        Primer::diff_3_end_flag = diff_3_end_flag;
    }

    const vector<DifferenceSite> &getDiff_vector() const {
        return diff_vector;
    }

    void setDiff_vector(const vector<DifferenceSite> &diff_vector) {
        Primer::diff_vector = diff_vector;
    }

    double getScore() const {
        return score;
    }

    void setScore(double score) {
        Primer::score = score;
    }

    int getStart_in_alignment() const {
        return start_in_alignment;
    }

    void setStart_in_alignment(int start_in_alignment) {
        Primer::start_in_alignment = start_in_alignment;
    }

    int getEnd_in_alignment() const {
        return end_in_alignment;
    }

    void setEnd_in_alignment(int end_in_alignment) {
        Primer::end_in_alignment = end_in_alignment;
    }

    int getDeletion_site_number() const {
        return deletion_site_number;
    }

    void setDeletion_site_number(int deletion_site_number) {
        Primer::deletion_site_number = deletion_site_number;
    }
};

class Primer_Pair {
private:
    Primer left_primer;
    Primer right_primer;
    double score;

public:
    const Primer &getLeft_primer() const {
        return left_primer;
    }

    void setLeft_primer(const Primer &left_primer) {
        Primer_Pair::left_primer = left_primer;
    }

    const Primer &getRight_primer() const {
        return right_primer;
    }

    void setRight_primer(const Primer &right_primer) {
        Primer_Pair::right_primer = right_primer;
    }

    double getScore() const {
        return score;
    }

    void setScore(double score) {
        Primer_Pair::score = score;
    }
};

bool fileExist(string path);
bool binExist(string path, string type);
bool flagExist(int argc, char** argv, string flag);
string generateString(int length);
string generateIDsFileFromMultipleFasta(string fasta_file_path);
void mergeAllPrimersResult(string query_id, string primers_result_path, string output_path);
map<string, string> generateIDsFileFromBlastTab(string blast_table, int id_hits_for_each_query);
string designPrimer(string id_file, string blast_table, string dbs, string bed_path, string muscle_path, string primer3_path, int flanking_length, int product_min_length, int product_max_length, int diff_setting, string end_flag);
map<string, string> readIDMap(string id_file);
map<string, int> readAllSeqLength(string fai_file_path);
string generateBedLocationFile(map<string, string> id_map, string blast_tab, string id_file, int flanking_length, string dbs);
void parseAlignment(string alignment_file_path);
int getMin(int n1, int n2, int n3, int n4);
int getMax(int n1, int n2, int n3, int n4);
std::vector<std::string> split(std::string str,std::string pattern);
string& replace_all(string&   str,const   string&   old_value,const   string&   new_value);
string generatePrimer3Input(string alignment_file_path, map<string, string> parameters);
map<string, string> getSeqMap(string fasta_file_path, string type);
map<string, vector<DifferenceSite>> identifyDifferentSite(map<string, string> alignment_seq_map);
vector<Primer> getSpecificPrimer(string primer3_output_path, map<string, string> seq_map, map<string, string> alignment_seq_map, map<string, vector<DifferenceSite>> ds_map, string type);
vector<Primer> parsePrimerVector(vector<Primer> primer_vector);
bool cmpScore(Primer p1, Primer p2);
bool cmpPairScore(Primer_Pair p1, Primer_Pair p2);
void generatePrimerPairResult(vector<Primer> primer_vector_left, vector<Primer> primer_vector_right, int min_product_size, int max_product_size, string filter_type, int total_diff, int terminal_diff, int terminal_len, string last_base_diff_flag);
void generatePrimerPairResult_old(vector<Primer> primer_vector_left, vector<Primer> primer_vector_right, int min_product_size, int max_product_size, int diff_setting, string end_flag, string output, string output_download);
int getAlignmentSite(string alignment_seq, int real_site);
string generateDifferentSiteStringInPrimer(Primer p);
string getReverseComplement(string seq);
int calculateSeqSitebyAlignmentSite(string seq, int alignment_site);
int calculateDiffSiteIn3Terminal(Primer p, int len, string type);
void calculatePrimerDimer(string primer3_output);

string input_multiple_fasta_path = "";
string parameters_file = "";

int main(int argc, char** argv) {

    if(argc < 2) {
        cout << "Usage: GSP" << "\n";
        //cout << "-i    id file path" << "\n";
        cout << "-r    blast table result path" << "\n";
        cout << "-d    sequence database path (blast database path)" << "\n";
        cout << "-a    fasta file path (for design specefic priemrs in multiple sequences only)" << "\n";
        cout << "-b    bedtools path (default: bedtools)" << "\n";
        cout << "-m    muscle path (default: muscle)" << "\n";
        cout << "-p    primer3 path (default: primer3_core)" << "\n";
        cout << "-t    primer3 parameters file path" << "\n";
        cout << "-o    output path" << "\n";
        cout << "-q    number of hit sequences for perimer design of each query (default: 3)" << "\n";
        cout << "-f    hit sequence flanking length (default: 200)" << "\n";
        cout << "-s    product min size (default: 200)" << "\n";
        cout << "-l    product max size (default: 1000)" << "\n";
        cout << "-c    different site in primer (default: 2)" << "\n";
        cout << "-e    different site in 3 end of primer (default: No)" << "\n";

        cout << "Note:" << "\n";
        cout << "1. Make sure that the dependencies (bedtools, muscle, primer3_core) are in the $PATH. If not, they must be specified (-b, -m, -p)." << "\n";
        cout << "2. Please use absolute path." << "\n";
        exit(0);
    }

    string id_file;
    string blast_table;
    string dbs = "";
    string bed_path = "bedtools";
    string muscle_path = "muscle";
    string primer3_path = "primer3_core";
    int id_hits_for_each_query = 3;
    int flanking_length = 200;
    int product_min_length = 200;
    int product_max_length = 1000;
    int diff_setting = 2;
    string end_flag = "No";
    string output_path = blast_table.substr(0, blast_table.find_last_of("/")) + "/" + "primers.csv";
    string web_flag = "No";

    int bs = 0;
    int ms = 0;
    int ps = 0;

    if(!flagExist(argc, argv, "-a")) {
        if(!flagExist(argc, argv, "-r") || !flagExist(argc, argv, "-d")) {
            if(!flagExist(argc, argv, "-r")) {
                cout << "Input error: (-r) blast table result path must be specified" << "\n";
                exit(0);
            }
            if(!flagExist(argc, argv, "-d")) {
                cout << "Input error: (-d) sequence database path (blast database path) must be specified" << "\n";
                exit(0);
            }
        }
    }

    if(!flagExist(argc, argv, "-o")) {
        cout << "Input error: (-o) output path must be specified" << "\n";
        exit(0);
    }

    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], "-i") == 0) {
            id_file = argv[i + 1];
        }
        if(strcmp(argv[i], "-r") == 0) {
            blast_table = argv[i + 1];
            if(!fileExist(blast_table)) {
                cout << "Input error: (-r) blast table result not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-d") == 0) {
            dbs = argv[i + 1];
            if(!fileExist(dbs)) {
                cout << "Input error: (-d) sequence database (blast database) not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-a") == 0) {
            input_multiple_fasta_path = argv[i + 1];
            if(!fileExist(input_multiple_fasta_path)) {
                cout << "Input error: (-a) fasta file not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-b") == 0) {
            bed_path = argv[i + 1];
            bs = 1;
            if(!binExist(bed_path, "s")) {
                cout << "Input error: (-b) bedtools not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-m") == 0) {
            muscle_path = argv[i + 1];
            ms = 1;
            if(!binExist(muscle_path, "s")) {
                cout << "Input error: (-m) muscle not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-p") == 0) {
            primer3_path = argv[i + 1];
            ps = 1;
            if(!binExist(primer3_path, "s")) {
                cout << "Input error: (-p) primer3_core not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-t") == 0) {
            parameters_file = argv[i + 1];
            if(!fileExist(parameters_file)) {
                cout << "Input error: (-t) primer3 parameters file not found" << "\n";
                exit(0);
            }
        }
        if(strcmp(argv[i], "-q") == 0) {
            id_hits_for_each_query = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-f") == 0) {
            flanking_length = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-s") == 0) {
            product_min_length = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-l") == 0) {
            product_max_length = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-c") == 0) {
            diff_setting = atoi(argv[i + 1]);
        }
        if(strcmp(argv[i], "-e") == 0) {
            end_flag = argv[i + 1];
        }
        if(strcmp(argv[i], "-o") == 0) {
            output_path = argv[i + 1];
        }
        if(strcmp(argv[i], "-w") == 0) {
            web_flag = argv[i + 1];
        }
    }

    if(bs == 0) {
        if(!binExist(bed_path, "")) {
            cout << "Input error: (-b) bedtools not found" << "\n";
            exit(0);
        }
    }
    if(ms == 0) {
        if(!binExist(muscle_path, "")) {
            cout << "Input error: (-m) muscle not found" << "\n";
            exit(0);
        }
    }
    if(ps == 0) {
        if(!binExist(primer3_path, "")) {
            cout << "Input error: (-p) primer3_core not found" << "\n";
            exit(0);
        }
    }

/*
    string blast_table = "/home/yiwang/Public/GSP/example/blast.tab";
    string dbs = "/var/www/GSP/dbs/Triticum_aestivum.fa";
    string bed_path = "bedtools";
    string muscle_path = "/var/www/GSP/bin/muscle3.8.31_i86linux64";
    string primer3_path = "/var/www/GSP/bin/primer3_core";
    int flanking_length = 200;
    int product_min_length = 200;
    int product_max_length = 1000;
    int diff_setting = 2;
    string end_flag = "No";

    string id_file = "/home/yiwang/Public/GSP/ids";
    string output_path = blast_table.substr(0, blast_table.find_last_of("/")) + "/" + "primers.csv";
    string web_flag = "No";
*/
/*
    string blast_table = "/home/yiwang/Public/GSP/example/blast.tab";
    string dbs = "/var/www/GSP/dbs/Triticum_aestivum.fa";
    string bed_path = "bedtools";
    string muscle_path = "/var/www/GSP/bin/muscle3.8.31_i86linux64";
    string primer3_path = "/var/www/GSP/bin/primer3_core";
    int flanking_length = 200;
    int product_min_length = 200;
    int product_max_length = 1000;
    int diff_setting = 2;
    string end_flag = "No";

    string id_file = "/home/yiwang/Public/GSP/ids";
    string output_path = blast_table.substr(0, blast_table.find_last_of("/")) + "/" + "primers.csv";
    string web_flag = "No";
*/
/*
    string blast_table = "/home/yiwang/Public/GSP/example/blast.tab";
    string dbs = "/var/www/GSP/dbs/Triticum_aestivum.fa";

    string web_flag = "Yes";
*/

    if(web_flag.compare("Yes") == 0) {
        designPrimer(id_file, blast_table, dbs, bed_path, muscle_path, primer3_path, flanking_length, product_min_length, product_max_length, diff_setting, end_flag);
    } else {
        if(input_multiple_fasta_path.compare("") == 0) {
            ifstream fin_exist(output_path);
            if(fin_exist) {
                remove(output_path.c_str());
            }
            fin_exist.close();

            cout << "Parse blast table result..." << "\n";

            map<string, string> ids_path_map = generateIDsFileFromBlastTab(blast_table, id_hits_for_each_query);

            cout << "Start design primer..." << "\n";

            map<string,string>::iterator iter;
            for(iter = ids_path_map.begin(); iter != ids_path_map.end(); ++iter) {
                string primers_result_path = designPrimer(iter->second, blast_table, dbs, bed_path, muscle_path, primer3_path, flanking_length, product_min_length, product_max_length, diff_setting, end_flag);
                mergeAllPrimersResult(iter->first, primers_result_path, output_path);
            }
        } else {
            ifstream fin_exist(output_path);
            if(fin_exist) {
                remove(output_path.c_str());
            }
            fin_exist.close();

            string ids_path = generateIDsFileFromMultipleFasta(input_multiple_fasta_path);

            cout << "Start design primer..." << "\n";

            string primers_result_path = designPrimer(ids_path, blast_table, dbs, bed_path, muscle_path, primer3_path, flanking_length, product_min_length, product_max_length, diff_setting, end_flag);
            mergeAllPrimersResult("Specific primers for " + input_multiple_fasta_path, primers_result_path, output_path);
        }
    }

    return 0;
}

bool fileExist(string path) {
    ifstream fin(path);
    int flag = 0;
    if(fin) {
        flag = 1;
    } else {
        flag = 0;
    }
    fin.close();
    if(flag == 1) {
        return true;
    } else {
        return false;
    }
}

bool binExist(string path, string type) {
    //cout << path << "\n";
    if(type == "s") {
        int flag = 0;
        ifstream fin(path);
        if (fin) {
            flag = 1;
        } else {
            flag = 0;
        }
        fin.close();
        if(flag == 1) {
            return true;
        } else {
            return false;
        }
    } else {
        char* env_p = getenv("PATH");
        string s_env_p = env_p;
        vector<string> str_v = split(s_env_p, ":");
        int flag = 0;
        for(int i = 0; i < str_v.size(); i++) {
            ifstream fin(str_v[i] + "/" + path);
            if (fin) {
                flag = 1;
                break;
            }
            fin.close();
        }
        if(flag == 1) {
            return true;
        } else {
            return false;
        }
    }

}

bool flagExist(int argc, char** argv, string flag) {
    int exist_flag = 0;
    for(int i = 1; i < argc; i++) {
        if(strcmp(argv[i], flag.c_str()) == 0) {
            exist_flag = 1;
            break;
        }
    }
    if(exist_flag == 1) {
        return true;
    } else {
        return false;
    }
}

string generateString(int length) {
    static const char alphabet[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    int strLen = sizeof(alphabet)-1;
    char genRandom;
    string str;
    srand(time(0));
    for(int i=0; i<length; i++)
    {
        genRandom = alphabet[rand()%strLen];
        str += genRandom;
    }
    return str;
}

string generateIDsFileFromMultipleFasta(string fasta_file_path) {
    string temp_path = fasta_file_path.substr(0, fasta_file_path.find_last_of("/")) + "/" + generateString(6);
    string mkdir_tmp = "mkdir " + temp_path;
    system(mkdir_tmp.c_str());

    string folder_path = temp_path + "/" + "1";
    string mkdir_cmd = "mkdir " + folder_path;
    system(mkdir_cmd.c_str());

    ofstream fout(folder_path + "/" + "ids");
    ifstream fin(fasta_file_path);
    string line;
    while(getline(fin, line)) {
        if(line.find(">") <= line.length()) {
            vector<string> tmp_v = split(line, "|");
            fout << replace_all(tmp_v[0], ">", "") << "\n";
        }
    }
    fin.close();
    fout.close();

    return folder_path + "/" + "ids";
}

map<string, string> readParametersFile(string parameters_file_path) {
    map<string, string> parameters;

    parameters["PRIMER_OPT_SIZE"]="20";
    parameters["PRIMER_MIN_SIZE"]="18";
    parameters["PRIMER_MAX_SIZE"]="23";
    parameters["PRIMER_PRODUCT_OPT_SIZE"]="200";
    parameters["PRIMER_PRODUCT_SIZE_RANGE"]="150-700";
    parameters["PRIMER_MIN_GC"]="30.0";
    parameters["PRIMER_MAX_GC"]="70.0";
    parameters["PRIMER_MIN_TM"]="57.0";
    parameters["PRIMER_OPT_TM"]="60.0";
    parameters["PRIMER_MAX_TM"]="63.0";
    parameters["PRIMER_MAX_END_STABILITY"]="20.0";
    parameters["PRIMER_PAIR_MAX_DIFF_TM"]="3.0";
    parameters["PRIMER_MAX_TEMPLATE_MISPRIMING"]="12";
    parameters["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"]="24";
    parameters["PRIMER_PAIR_MAX_COMPL_ANY"]="1.00";

    if(parameters_file_path.compare("") != 0) {
        ifstream fin(parameters_file_path);
        string line;
        while(getline(fin, line)) {
            if(line.find("PRIMER_OPT_SIZE") <= line.length()) {
                parameters["PRIMER_OPT_SIZE"] = replace_all(line, "PRIMER_OPT_SIZE=", "");
            }
            if(line.find("PRIMER_MIN_SIZE") <= line.length()) {
                parameters["PRIMER_MIN_SIZE"] = replace_all(line, "PRIMER_MIN_SIZE=", "");
            }
            if(line.find("PRIMER_MAX_SIZE") <= line.length()) {
                parameters["PRIMER_MAX_SIZE"] = replace_all(line, "PRIMER_MAX_SIZE=", "");
            }
            if(line.find("PRIMER_PRODUCT_OPT_SIZE") <= line.length()) {
                parameters["PRIMER_PRODUCT_OPT_SIZE"] = replace_all(line, "PRIMER_PRODUCT_OPT_SIZE=", "");
            }
            if(line.find("PRIMER_PRODUCT_SIZE_RANGE") <= line.length()) {
                parameters["PRIMER_PRODUCT_SIZE_RANGE"] = replace_all(line, "PRIMER_PRODUCT_SIZE_RANGE=", "");
            }
            if(line.find("PRIMER_MIN_GC") <= line.length()) {
                parameters["PRIMER_MIN_GC"] = replace_all(line, "PRIMER_MIN_GC=", "");
            }
            if(line.find("PRIMER_MAX_GC") <= line.length()) {
                parameters["PRIMER_MAX_GC"] = replace_all(line, "PRIMER_MAX_GC=", "");
            }
            if(line.find("PRIMER_MIN_TM") <= line.length()) {
                parameters["PRIMER_MIN_TM"] = replace_all(line, "PRIMER_MIN_TM=", "");
            }
            if(line.find("PRIMER_OPT_TM") <= line.length()) {
                parameters["PRIMER_OPT_TM"] = replace_all(line, "PRIMER_OPT_TM=", "");
            }
            if(line.find("PRIMER_MAX_TM") <= line.length()) {
                parameters["PRIMER_MAX_TM"] = replace_all(line, "PRIMER_MAX_TM=", "");
            }
            if(line.find("PRIMER_MAX_END_STABILITY") <= line.length()) {
                parameters["PRIMER_MAX_END_STABILITY"] = replace_all(line, "PRIMER_MAX_END_STABILITY=", "");
            }
            if(line.find("PRIMER_PAIR_MAX_DIFF_TM") <= line.length()) {
                parameters["PRIMER_PAIR_MAX_DIFF_TM"] = replace_all(line, "PRIMER_PAIR_MAX_DIFF_TM=", "");
            }
            if(line.find("PRIMER_MAX_TEMPLATE_MISPRIMING") <= line.length()) {
                parameters["PRIMER_MAX_TEMPLATE_MISPRIMING"] = replace_all(line, "PRIMER_MAX_TEMPLATE_MISPRIMING=", "");
            }
            if(line.find("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING") <= line.length()) {
                parameters["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"] = replace_all(line, "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=", "");
            }
            if(line.find("PRIMER_PAIR_MAX_COMPL_ANY") <= line.length()) {
                parameters["PRIMER_PAIR_MAX_COMPL_ANY"] = replace_all(line, "PRIMER_PAIR_MAX_COMPL_ANY=", "");
            }
        }
        fin.close();
    }

    return  parameters;
};

void mergeAllPrimersResult(string query_id, string primers_result_path, string output_path) {
    ofstream fout(output_path, ios::app);
    fout << "# Specific primers for " + query_id << "\n";
    ifstream fin(primers_result_path);
    string line;
    while(getline(fin, line)) {
        fout << line << "\n";
    }
    fin.close();
    fout << "\n";
    fout.close();
}

map<string, string> generateIDsFileFromBlastTab(string blast_table, int id_hits_for_each_query) {
    string temp_path = blast_table.substr(0, blast_table.find_last_of("/")) + "/" + generateString(6);
    string mkdir_tmp = "mkdir " + temp_path;
    system(mkdir_tmp.c_str());

    //vector<string> ids_file_vec;
    map<string, string> ids_map;

    ifstream fin(blast_table);
    string line;
    string query_id = "";
    string hits_id = "";
    int hit_number = 0;
    int index = 0;
    map<string, string> hits_map;
    while(getline(fin, line)) {
        vector<string> tmp_v = split(line, "\t");
        if(tmp_v[0].compare(query_id) != 0) {
            if(query_id.compare("") != 0) {
                if(hit_number == id_hits_for_each_query) {
                    index++;
                    stringstream ss;
                    ss << index;
                    string folder_path = temp_path + "/" + ss.str().c_str();
                    string mkdir_cmd = "mkdir " + folder_path;
                    system(mkdir_cmd.c_str());

                    ofstream fout(folder_path + "/" + "ids");
                    fout << hits_id;
                    fout.close();
                    //ids_file_vec.push_back(folder_path + "/" + "ids");
                    ids_map[query_id] = folder_path + "/" + "ids";
                }
            }
            query_id = tmp_v[0];
            hits_id = "";
            hits_id = hits_id + tmp_v[1] + "\n";
            hit_number = 1;
            hits_map.clear();
            hits_map[tmp_v[1]] = "";
        } else {
            if(hits_map.find(tmp_v[1]) == hits_map.end()) {
                hits_map[tmp_v[1]] = "";
                hit_number++;
                if(hit_number > id_hits_for_each_query) {
                    hit_number = id_hits_for_each_query;
                    continue;
                }
                hits_id = hits_id + tmp_v[1] + "\n";
            }
        }
    }

    if(hit_number == id_hits_for_each_query) {
        index++;
        stringstream ss;
        ss << index;
        string folder_path = temp_path + "/" + ss.str().c_str();
        string mkdir_cmd = "mkdir " + folder_path;
        system(mkdir_cmd.c_str());

        ofstream fout(folder_path + "/" + "ids");
        fout << hits_id;
        fout.close();
        //ids_file_vec.push_back(folder_path + "/" + "ids");
        ids_map[query_id] = folder_path + "/" + "ids";
    }

    fin.close();
    return ids_map;
}

string designPrimer(string id_file, string blast_table, string dbs, string bed_path, string muscle_path, string primer3_path, int flanking_length, int product_min_length, int product_max_length, int diff_setting, string end_flag) {
    string alignment_file_path = id_file.substr(0, id_file.find_last_of("/")) + "/alignment.fas";
    map<string, string> id_map = readIDMap(id_file);
    if(dbs.compare("") != 0) {
        string bed_file = generateBedLocationFile(id_map, blast_table, id_file, flanking_length, dbs);
        string extract_fas = id_file.substr(0, id_file.find_last_of("/")) + "/extract.fas";
        string bed_cmd = bed_path + " getfasta -fi " + dbs + " -bed " + bed_file + " -s -fo " + extract_fas;
        cout << bed_cmd << "\n";
        system(bed_cmd.c_str());
        string align_cmd = muscle_path + " -in " + extract_fas + " -out " + alignment_file_path + " -quiet";
        system(align_cmd.c_str());
    } else {
        string tmp_str = id_file.substr(0, id_file.find_last_of("/"));

        string all_fas;
        if(input_multiple_fasta_path.compare("") == 0) {
            all_fas = tmp_str.substr(0, tmp_str.find_last_of("/")) + "/msa_input.seq";
        } else {
            all_fas = input_multiple_fasta_path;
        }

        string extract_fas = tmp_str + "/extract.fas";

        ifstream fin(all_fas);
        ofstream fout(tmp_str + "/extract.fas");
        string line;
        string id = "";
        string seq;
        while (getline(fin, line)) {
            boost::algorithm::trim(line);
            if (line.find(">") < line.length()) {
                if (id.compare("") != 0) {
                    if (id_map.find(id) != id_map.end()) {
                        fout << ">" + id + "\n" + seq << "\n";
                    }
                }
                id = replace_all(line, ">", "");
                seq = "";
            } else {
                seq += line;
            }
        }
        if (id_map.find(id) != id_map.end()) {
            fout << ">" + id + "\n" + seq << "\n";
        }
        fin.close();
        fout.close();

        string align_cmd = muscle_path + " -in " + extract_fas + " -out " + alignment_file_path + " -quiet";
        system(align_cmd.c_str());
    }
    parseAlignment(alignment_file_path);

    string alignment_file_format_path = id_file.substr(0, id_file.find_last_of("/")) + "/alignment_format.fas";

    map<string, string> parameters = readParametersFile(parameters_file);
    string primer3_input_path = generatePrimer3Input(alignment_file_format_path, parameters);

    ifstream primer3_output_fin(primer3_input_path.substr(0, primer3_input_path.find_last_of("/")) + "/" + "primer3.output");
    string primer3_cmd = primer3_path + " -default_version=1 -format_output -output=" + primer3_input_path.substr(0, primer3_input_path.find_last_of("/")) + "/" + "primer3.output" + " " + primer3_input_path;
    system(primer3_cmd.c_str());

    //calculatePrimerDimer(primer3_input_path.substr(0, primer3_input_path.find_last_of("/")) + "/" + "primer3.output");

    map<string, string> seq_map = getSeqMap(alignment_file_format_path, "alignment"); //remove dash in alignment, get normal sequence
    map<string, string> alignment_seq_map = getSeqMap(alignment_file_format_path, "normal"); //get alignment sequence
    map<string, vector<DifferenceSite>> ds_map = identifyDifferentSite(alignment_seq_map);

    vector<Primer> primer_vector_left = getSpecificPrimer(primer3_input_path.substr(0, primer3_input_path.find_last_of("/")) + "/" + "primer3.output", seq_map, alignment_seq_map, ds_map, "LEFT_PRIMER");
    vector<Primer> primer_vector_right = getSpecificPrimer(primer3_input_path.substr(0, primer3_input_path.find_last_of("/")) + "/" + "primer3.output", seq_map, alignment_seq_map, ds_map, "RIGHT_PRIMER");

    //cout << primer_vector_left.size() << "\n";
    //cout << primer_vector_right.size() << "\n";

    vector<Primer> primer_vector_left_p = parsePrimerVector(primer_vector_left);
    vector<Primer> primer_vector_right_p = parsePrimerVector(primer_vector_right);

    //generatePrimerPairResult(primer_vector_left, primer_vector_right, 150, 1000, "ABB", 1, 0, 3, "Yes");
    generatePrimerPairResult_old(primer_vector_left_p, primer_vector_right_p, product_min_length, product_max_length, diff_setting, end_flag, alignment_file_format_path.substr(0, alignment_file_format_path.find_last_of("/")) + "/primer_pairs.tab", alignment_file_format_path.substr(0, alignment_file_format_path.find_last_of("/")) + "/all_primer_pairs.csv");

    return alignment_file_format_path.substr(0, alignment_file_format_path.find_last_of("/")) + "/all_primer_pairs.csv";
}

map<string, string> readIDMap(string id_file) {
    map<string, string> id_map;
    ifstream fin(id_file);
    string line;
    while(getline(fin, line)) {
        //cout << line << endl;
        id_map.insert(pair<string, string>(line, ""));
    }
    fin.close();
    return id_map;
}

std::vector<std::string> split(std::string str,std::string pattern)
{
    std::string::size_type pos;
    std::vector<std::string> result;
    str+=pattern;//扩展字符串以方便操作
    int size=str.size();

    for(int i=0; i<size; i++)
    {
        pos=str.find(pattern,i);
        if(pos<size)
        {
            std::string s=str.substr(i,pos-i);
            result.push_back(s);
            i=pos+pattern.size()-1;
        }
    }
    return result;
}

map<string, int> readAllSeqLength(string fai_file_path) {
    map<string, int> seq_length_map;
    ifstream fin(fai_file_path);
    string line;
    while(getline(fin, line)) {
        vector<string> tmp_v = split(line, "\t");
        seq_length_map[tmp_v[0]] = atoi(tmp_v[1].c_str());
    }
    fin.close();
    return seq_length_map;
};

string generateBedLocationFile(map<string, string> id_map, string blast_tab, string id_file, int flanking_length, string dbs) {
    ifstream fin(blast_tab);
    int pos = id_file.find_last_of("/");
    ofstream fout(id_file.substr(0, pos) + "/bed");
    string line;
    string id = "";
    int first_hit_start = 0;
    int first_hit_end = 0;
    int min_hit_start = 0;
    int max_hit_start = 0;
    int min_hit_end = 0;
    int max_hit_end = 0;
    string strand = "";

    map<string, int> seq_length_map = readAllSeqLength(dbs + ".fai");

    while(getline(fin, line)) {
        vector<string> str_vector;
        boost::split(str_vector, line, boost::is_any_of("\t"), boost::token_compress_on);
        map<string, string>::iterator it;
        it = id_map.find(str_vector[1]);
        if(it != id_map.end()) {
            if (str_vector[1].compare(id) != 0) {
                if (id.compare("") != 0) {
                    int tmp_min = getMin(min_hit_start, min_hit_end, max_hit_start, max_hit_end) - flanking_length;
                    int tmp_max = getMax(min_hit_start, min_hit_end, max_hit_start, max_hit_end) + flanking_length;
                    if(tmp_max - tmp_min > 30000) {
                        if(first_hit_start < first_hit_end) {
                            int tmp_first_start = first_hit_start - flanking_length;
                            if(tmp_first_start <= 0) {
                                tmp_first_start = 1;
                            }
                            int tmp_first_end = first_hit_end + flanking_length;
                            if(seq_length_map.find(id) != seq_length_map.end()) {
                                int len = seq_length_map[id];
                                if(tmp_first_end > len) {
                                    tmp_first_end = len;
                                }
                            }
                            fout << id << "\t" << tmp_first_start << "\t" << tmp_first_end << "\t" << id << "\t" << "1" << "\t" << strand << "\n";
                        } else {
                            int tmp_first_start = first_hit_end - flanking_length;
                            if(tmp_first_start <= 0) {
                                tmp_first_start = 1;
                            }
                            int tmp_first_end = first_hit_start + flanking_length;
                            if(seq_length_map.find(id) != seq_length_map.end()) {
                                int len = seq_length_map[id];
                                if(tmp_first_end > len) {
                                    tmp_first_end = len;
                                }
                            }
                            fout << id << "\t" << tmp_first_start << "\t" << tmp_first_end << "\t" << id << "\t" << "1" << "\t" << strand << "\n";
                        }
                    } else {
                        if(tmp_min <= 0) {
                            tmp_min = 1;
                        }
                        if(seq_length_map.find(id) != seq_length_map.end()) {
                            int len = seq_length_map[id];
                            if(tmp_max > len) {
                                tmp_max = len;
                            }
                        }
                        fout << id << "\t" << tmp_min << "\t" << tmp_max << "\t" << id << "\t" << "1" << "\t" << strand << "\n";
                    }
                }
                id = str_vector[1];
                first_hit_start = atoi(str_vector[8].c_str());
                first_hit_end = atoi(str_vector[9].c_str());
                min_hit_start = atoi(str_vector[8].c_str());
                max_hit_start = atoi(str_vector[8].c_str());
                min_hit_end = atoi(str_vector[9].c_str());
                max_hit_end = atoi(str_vector[9].c_str());
                if(min_hit_start < min_hit_end) {
                    strand = "+";
                } else {
                    strand = "-";
                }
            } else {
                if (min_hit_start - atoi(str_vector[8].c_str()) <= 5000 && min_hit_start - atoi(str_vector[8].c_str()) > 0) {
                    min_hit_start = atoi(str_vector[8].c_str());
                }
                if (atoi(str_vector[8].c_str()) - max_hit_start <= 5000 && atoi(str_vector[8].c_str()) - max_hit_start > 0) {
                    max_hit_start = atoi(str_vector[8].c_str());
                }
                if (min_hit_end - atoi(str_vector[9].c_str()) <= 5000 && min_hit_end - atoi(str_vector[9].c_str()) > 0) {
                    min_hit_end = atoi(str_vector[9].c_str());
                }
                if (atoi(str_vector[9].c_str()) - max_hit_end <= 5000 && atoi(str_vector[9].c_str()) - max_hit_end > 0) {
                    max_hit_end = atoi(str_vector[9].c_str());
                }
            }
        }
    }

    if (id.compare("") != 0) {
        int tmp_min = getMin(min_hit_start, min_hit_end, max_hit_start, max_hit_end) - flanking_length;
        int tmp_max = getMax(min_hit_start, min_hit_end, max_hit_start, max_hit_end) + flanking_length;
        if(tmp_max - tmp_min > 30000) {
            if(first_hit_start < first_hit_end) {
                int tmp_first_start = first_hit_start - flanking_length;
                if(tmp_first_start <= 0) {
                    tmp_first_start = 1;
                }
                int tmp_first_end = first_hit_end + flanking_length;
                if(seq_length_map.find(id) != seq_length_map.end()) {
                    int len = seq_length_map[id];
                    if(tmp_first_end > len) {
                        tmp_first_end = len;
                    }
                }
                fout << id << "\t" << tmp_first_start << "\t" << tmp_first_end << "\t" << id << "\t" << "1" << "\t" << strand << "\n";
            } else {
                int tmp_first_start = first_hit_end - flanking_length;
                if(tmp_first_start <= 0) {
                    tmp_first_start = 1;
                }
                int tmp_first_end = first_hit_start + flanking_length;
                if(seq_length_map.find(id) != seq_length_map.end()) {
                    int len = seq_length_map[id];
                    if(tmp_first_end > len) {
                        tmp_first_end = len;
                    }
                }
                fout << id << "\t" << tmp_first_start << "\t" << tmp_first_end << "\t" << id << "\t" << "1" << "\t" << strand << "\n";
            }
        } else {
            if(tmp_min <= 0) {
                tmp_min = 1;
            }
            if(seq_length_map.find(id) != seq_length_map.end()) {
                int len = seq_length_map[id];
                if(tmp_max > len) {
                    tmp_max = len;
                }
            }
            fout << id << "\t" << tmp_min << "\t" << tmp_max << "\t" << id << "\t" << "1" << "\t" << strand << "\n";
        }
    }

    fin.close();
    fout.close();

    return id_file.substr(0, pos) + "/bed";
}

int getMin(int n1, int n2, int n3, int n4) {
    int a,b;
    if(n1 < n2) {
        a = n1;
    } else {
        a = n2;
    }

    if(n3 < n4) {
        b = n3;
    } else {
        b = n4;
    }

    if(a < b) {
        return a;
    } else {
        return b;
    }
}

int getMax(int n1, int n2, int n3, int n4) {
    int a,b;
    if(n1 > n2) {
        a = n1;
    } else {
        a = n2;
    }

    if(n3 > n4) {
        b = n3;
    } else {
        b = n4;
    }

    if(a > b) {
        return a;
    } else {
        return b;
    }
}

void parseAlignment(string alignment_file_path) {
    string output_path = alignment_file_path.substr(0, alignment_file_path.find_last_of("/")) + "/" + "alignment_format.fas";
    map<string, string> alignment_map = getSeqMap(alignment_file_path, "normal");
    map<string, string>::iterator iter = alignment_map.begin();
    int length = 0;
    for(; iter != alignment_map.end(); ++iter) {
        string seq = iter->second;
        //cout << seq << endl;
        length = seq.length();
        break;
    }

    vector<string> con_vector;
    for(int i = 0; i < length; i++) {
        string concensus;
        map<string, string>::iterator it = alignment_map.begin();
        for(; it != alignment_map.end(); ++it) {
            string seq = it->second;
            concensus += seq.at(i);
        }
        //cout << concensus << endl;
        con_vector.push_back(concensus);
    }

    int start_a = 0;
    for(int i = 0; i < con_vector.size(); i++) {
        string con = con_vector[i];
        int dash_number = 0;
        for(int j = 0; j < con.length(); j++) {
            if(con.at(j) == '-') {
                dash_number++;
            }
        }
        if(dash_number <= con.length() - 2) {
            start_a = i + 1;
            break;
        }
    }

    int end_a = 0;
    for(int i = con_vector.size() - 1; i >= 0; i--) {
        string con = con_vector[i];
        int dash_number = 0;
        for(int j = 0; j < con.length(); j++) {
            if(con.at(j) == '-') {
                dash_number++;
            }
        }
        if(dash_number <= con.length() - 2) {
            end_a = i + 1;
            break;
        }
    }

    //cout << start_a << " " << end_a << endl;

    ifstream fin(alignment_file_path);
    ofstream fout(output_path);
    string line;
    string id;
    string seq;
    while(getline(fin, line)) {
        if(line.find(">") < line.length()) {
            if(id.compare("") != 0) {
                fout << ">" + id + "\n" << seq.substr(start_a - 1, end_a - start_a + 1) << "\n";
            }
            id = replace_all(line, ">", "");
            seq = "";
        } else {
            seq += line;
        }
    }

    fout << ">" + id + "\n" << seq.substr(start_a - 1, end_a - start_a + 1) << "\n";

    fin.close();
    fout.close();
}

string generatePrimer3Input(string alignment_file_path, map<string, string> parameters) {
    ifstream fin(alignment_file_path);
    string primer3_input_path = alignment_file_path.substr(0, alignment_file_path.find_last_of("/")) + "/" + "primer3.input";
    ofstream fout(primer3_input_path);
    string line;
    string id;
    string aligned_seq;
    while(getline(fin, line)) {
        if(line.find(">") < line.length()) {
            if(id.compare("") != 0) {
                fout << "SEQUENCE_ID=" << id << "\n";
                fout << "SEQUENCE_TEMPLATE=" << replace_all(aligned_seq, "-", "") << "\n";
                fout << "PRIMER_TASK=" << "pick_primer_list" << "\n";
                fout << "PRIMER_MAX_NS_ACCEPTED=" << "0" << "\n";

                fout << "PRIMER_OPT_SIZE=" << parameters["PRIMER_OPT_SIZE"] << "\n";
                fout << "PRIMER_MIN_SIZE=" << parameters["PRIMER_MIN_SIZE"] << "\n";
                fout << "PRIMER_MAX_SIZE=" << parameters["PRIMER_MAX_SIZE"] << "\n";
                fout << "PRIMER_PRODUCT_OPT_SIZE=" << parameters["PRIMER_PRODUCT_OPT_SIZE"] << "\n";
                fout << "PRIMER_PRODUCT_SIZE_RANGE=" << parameters["PRIMER_PRODUCT_SIZE_RANGE"] << "\n";
                fout << "PRIMER_MIN_GC=" << parameters["PRIMER_MIN_GC"] << "\n";
                fout << "PRIMER_MAX_GC=" << parameters["PRIMER_MAX_GC"] << "\n";
                fout << "PRIMER_MIN_TM=" << parameters["PRIMER_MIN_TM"] << "\n";
                fout << "PRIMER_OPT_TM=" << parameters["PRIMER_OPT_TM"] << "\n";
                fout << "PRIMER_MAX_TM=" << parameters["PRIMER_MAX_TM"] << "\n";
                fout << "PRIMER_MAX_END_STABILITY=" << parameters["PRIMER_MAX_END_STABILITY"] << "\n";
                fout << "PRIMER_PAIR_MAX_DIFF_TM=" << parameters["PRIMER_PAIR_MAX_DIFF_TM"] << "\n";
                fout << "PRIMER_MAX_TEMPLATE_MISPRIMING=" << parameters["PRIMER_MAX_TEMPLATE_MISPRIMING"] << "\n";
                fout << "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=" << parameters["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"] << "\n";
                fout << "PRIMER_PAIR_MAX_COMPL_ANY=" << parameters["PRIMER_PAIR_MAX_COMPL_ANY"] << "\n";
                /*
                fout << "PRIMER_OPT_SIZE=" << "20" << "\n";
                fout << "PRIMER_MIN_SIZE=" << "18" << "\n";
                fout << "PRIMER_MAX_SIZE=" << "23" << "\n";
                fout << "PRIMER_PRODUCT_OPT_SIZE=" << "200" << "\n";
                fout << "PRIMER_PRODUCT_SIZE_RANGE=" << "150-700" << "\n";
                fout << "PRIMER_MIN_GC=" << "30.0" << "\n";
                fout << "PRIMER_MAX_GC=" << "70.0" << "\n";
                fout << "PRIMER_MIN_TM=" << "57.0" << "\n";
                fout << "PRIMER_OPT_TM=" << "60.0" << "\n";
                fout << "PRIMER_MAX_TM=" << "63.0" << "\n";
                fout << "PRIMER_MAX_END_STABILITY=" << "20.0" << "\n";
                fout << "PRIMER_PAIR_MAX_DIFF_TM=" << "3.0" << "\n";
                fout << "PRIMER_MAX_TEMPLATE_MISPRIMING=" << "12" << "\n";
                fout << "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=" << "24" << "\n";
                fout << "PRIMER_PAIR_MAX_COMPL_ANY=" << "1.00" << "\n";
                */
                fout << "PRIMER_NUM_RETURN=" << "10000" << "\n";
                fout << "=" << "\n";
            }
            id = line.substr(1, line.length());
            aligned_seq = "";
        } else {
            aligned_seq += line;
        }
    }

    if(id.compare("") != 0) {
        fout << "SEQUENCE_ID=" << id << "\n";
        fout << "SEQUENCE_TEMPLATE=" << replace_all(aligned_seq, "-", "") << "\n";
        fout << "PRIMER_TASK=" << "pick_primer_list" << "\n";
        fout << "PRIMER_MAX_NS_ACCEPTED=" << "0" << "\n";

        fout << "PRIMER_OPT_SIZE=" << parameters["PRIMER_OPT_SIZE"] << "\n";
        fout << "PRIMER_MIN_SIZE=" << parameters["PRIMER_MIN_SIZE"] << "\n";
        fout << "PRIMER_MAX_SIZE=" << parameters["PRIMER_MAX_SIZE"] << "\n";
        fout << "PRIMER_PRODUCT_OPT_SIZE=" << parameters["PRIMER_PRODUCT_OPT_SIZE"] << "\n";
        fout << "PRIMER_PRODUCT_SIZE_RANGE=" << parameters["PRIMER_PRODUCT_SIZE_RANGE"] << "\n";
        fout << "PRIMER_MIN_GC=" << parameters["PRIMER_MIN_GC"] << "\n";
        fout << "PRIMER_MAX_GC=" << parameters["PRIMER_MAX_GC"] << "\n";
        fout << "PRIMER_MIN_TM=" << parameters["PRIMER_MIN_TM"] << "\n";
        fout << "PRIMER_OPT_TM=" << parameters["PRIMER_OPT_TM"] << "\n";
        fout << "PRIMER_MAX_TM=" << parameters["PRIMER_MAX_TM"] << "\n";
        fout << "PRIMER_MAX_END_STABILITY=" << parameters["PRIMER_MAX_END_STABILITY"] << "\n";
        fout << "PRIMER_PAIR_MAX_DIFF_TM=" << parameters["PRIMER_PAIR_MAX_DIFF_TM"] << "\n";
        fout << "PRIMER_MAX_TEMPLATE_MISPRIMING=" << parameters["PRIMER_MAX_TEMPLATE_MISPRIMING"] << "\n";
        fout << "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=" << parameters["PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING"] << "\n";
        fout << "PRIMER_PAIR_MAX_COMPL_ANY=" << parameters["PRIMER_PAIR_MAX_COMPL_ANY"] << "\n";
        /*
        fout << "PRIMER_OPT_SIZE=" << "20" << "\n";
        fout << "PRIMER_MIN_SIZE=" << "18" << "\n";
        fout << "PRIMER_MAX_SIZE=" << "23" << "\n";
        fout << "PRIMER_PRODUCT_OPT_SIZE=" << "200" << "\n";
        fout << "PRIMER_PRODUCT_SIZE_RANGE=" << "150-700" << "\n";
        fout << "PRIMER_MIN_GC=" << "30.0" << "\n";
        fout << "PRIMER_MAX_GC=" << "70.0" << "\n";
        fout << "PRIMER_MIN_TM=" << "57.0" << "\n";
        fout << "PRIMER_OPT_TM=" << "60.0" << "\n";
        fout << "PRIMER_MAX_TM=" << "63.0" << "\n";
        fout << "PRIMER_MAX_END_STABILITY=" << "20.0" << "\n";
        fout << "PRIMER_PAIR_MAX_DIFF_TM=" << "3.0" << "\n";
        fout << "PRIMER_MAX_TEMPLATE_MISPRIMING=" << "12" << "\n";
        fout << "PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=" << "24" << "\n";
        fout << "PRIMER_PAIR_MAX_COMPL_ANY=" << "1.00" << "\n";
        */
        fout << "PRIMER_NUM_RETURN=" << "10000" << "\n";
        fout << "=" << "\n";
    }

    fin.close();
    fout.close();

    return primer3_input_path;
}

map<string, string> getSeqMap(string fasta_file_path, string type) {
    map<string, string> seq_map;
    ifstream fin(fasta_file_path);
    string line;
    string id;
    string seq;
    while(getline(fin, line)) {
        if(line.find(">") < line.length()) {
            if(id.compare("") != 0) {
                seq_map.insert(pair<string, string>(id, seq));
            }
            id = replace_all(line, ">", "");
            seq = "";
        } else {
            if(type.compare("normal") == 0) {
                seq += line;
            } else if(type.compare("alignment") == 0) {
                seq += replace_all(line, "-", "");
            }
        }
    }
    seq_map.insert(pair<string, string>(id, seq));
    fin.close();
    return seq_map;
}

bool cmpScore(Primer p1, Primer p2) {
    if(p1.getScore() >  p2.getScore()) {
        return true;
    } else {
        return false;
    }
}

bool cmpPairScore(Primer_Pair p1, Primer_Pair p2) {
    if(p1.getScore() >  p2.getScore()) {
        return true;
    } else {
        return false;
    }
}

vector<Primer> getSpecificPrimer(string primer3_output_path, map<string, string> seq_map, map<string, string> alignment_seq_map, map<string, vector<DifferenceSite>> ds_map, string type) {
    vector<Primer> primer_vector;
    ifstream fin(primer3_output_path);
    string line;
    int flag = 0;
    string seq_id;
    while(getline(fin, line)) {
        if(line.find("PRIMER PICKING RESULTS FOR ") < line.length()) {
            flag = 1;
            seq_id = replace_all(line, "PRIMER PICKING RESULTS FOR ", "");
        }
        if(flag == 1 && line.find(" " + type) < line.length()) {
            boost::algorithm::trim(line);
            boost::xpressive::sregex rgx = sregex::compile("\\s+");
            boost::xpressive::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
            boost::xpressive::sregex_token_iterator end;

            Primer p;
            p.setId(seq_id);
            int index = 0;
            while(iter != end) {
                //cout << *iter << "\n";
                if(index == 2) {
                    string start_str = *iter;
                    p.setStart(atoi(start_str.c_str()) + 1);
                    p.setStart_in_alignment(getAlignmentSite(alignment_seq_map[seq_id], p.getStart()));
                } else if(index == 3) {
                    string start_str = *iter;
                    p.setLength(atoi(start_str.c_str()));
                    if(type.compare("LEFT_PRIMER") == 0) {
                        p.setEnd_in_alignment(getAlignmentSite(alignment_seq_map[seq_id], p.getStart() + p.getLength() - 1));
                    } else {
                        p.setEnd_in_alignment(getAlignmentSite(alignment_seq_map[seq_id], p.getStart() - p.getLength() + 1));
                    }
                } else if(index == 4) {
                    p.setTm(*iter);
                } else if(index == 5) {
                    p.setGc(*iter);
                } else if(index == 6) {
                    p.setAny(*iter);
                } else if(index == 7) {
                    p.setThree(*iter);
                } else if(index == 8) {
                    p.setSeq(*iter);
                    break;
                }

                index++;
                iter++;
            }

            int primer_flag = 0;
            for(map<string,string>::iterator it = seq_map.begin(); it!= seq_map.end(); ++it) {
                if(seq_id.compare(it->first) != 0) {
                    string temp_seq = it->second;
                    if(type.compare("LEFT_PRIMER") == 0) {
                        if(temp_seq.find(*iter) < temp_seq.length()) {
                            primer_flag = 1;
                            break;
                        }
                    } else {
                        string rs_primer = getReverseComplement(*iter);
                        if(temp_seq.find(rs_primer) < temp_seq.length()) {
                            primer_flag = 1;
                            break;
                        }
                    }
                }
            }
            if(primer_flag == 0) {
                primer_vector.push_back(p);
            }
        }
    }
    fin.close();

    for(int r = 0; r < primer_vector.size(); r++) {
        Primer p = primer_vector[r];
        vector<DifferenceSite> ds_vector = ds_map[p.getId()];
        vector<DifferenceSite> diff_vector;
        double score = 0.0;
        //cout << p.getId() << " " << ds_vector.size() << "\n";
        int diff_count = 0;
        int deletion_count = 0;
        string diff_3_end_flag = "No";
        for(int s = 0; s < ds_vector.size(); s++) {
            if(type.compare("LEFT_PRIMER") == 0) {
                //cout << p.getStart() << " " << ds_vector[s].getIn_seq_site() << "\n";
                if(p.getStart() + p.getLength() - 1 < ds_vector[s].getIn_seq_site()) {
                    break;
                }
                if(p.getStart() + p.getLength() - 1 == ds_vector[s].getIn_seq_site()) {
                    if(ds_vector[s].getChar_in_different() == '-') {
                        break;
                    }
                }
                if(p.getStart() <= ds_vector[s].getIn_seq_site()) {
                    diff_vector.push_back(ds_vector[s]);
                    diff_count++;
                    if(ds_vector[s].getChar_in_different() == '-') {
                        deletion_count++;
                    }
                    score += 100.0 / (p.getStart() + p.getLength() - ds_vector[s].getIn_seq_site());
                }
                if(p.getStart() + p.getLength() - 1 == ds_vector[s].getIn_seq_site()) {
                    diff_3_end_flag = "Yes";
                }
            } else {
                if(p.getStart() < ds_vector[s].getIn_seq_site()) {
                    break;
                }
                if(p.getStart() == ds_vector[s].getIn_seq_site()) {
                    if(ds_vector[s].getChar_in_different() == '-') {
                        break;
                    }
                }
                if(p.getStart() - p.getLength() + 1 <= ds_vector[s].getIn_seq_site()) {
                    diff_vector.push_back(ds_vector[s]);
                    diff_count++;
                    if(ds_vector[s].getChar_in_different() == '-') {
                        deletion_count++;
                    }
                    score += 100.0 / (ds_vector[s].getIn_seq_site() - p.getStart() + p.getLength());
                }
                if(p.getStart() - p.getLength() + 1 == ds_vector[s].getIn_seq_site()) {
                    diff_3_end_flag = "Yes";
                }
            }
        }
        p.setDiff_site_number(diff_count);
        p.setDeletion_site_number(deletion_count);
        p.setDiff_3_end_flag(diff_3_end_flag);
        p.setDiff_vector(diff_vector);
        //cout << p.getDiff_site_number() << "\n";
        p.setScore(score);
        primer_vector[r] = p;
    }

    sort(primer_vector.begin(), primer_vector.end(), cmpScore);
    return primer_vector;
}

vector<Primer> parsePrimerVector(vector<Primer> primer_vector) {
    map<string, vector<int>> exist_sites;
    vector<Primer> primer_vector_new;
    for(int i = 0; i < primer_vector.size(); i++) {
        Primer p = primer_vector[i];

        int exist_flag = 0;
        if(exist_sites.find(p.getId()) != exist_sites.end()) {
            vector<int> sites_vector = exist_sites[p.getId()];
            for(int j = 0; j < sites_vector.size(); j++) {
                if((p.getStart() - sites_vector[j] <= 10 && p.getStart() - sites_vector[j] >= 0) || (sites_vector[j] - p.getStart() <= 10 && sites_vector[j] - p.getStart() >= 0)) {
                    exist_flag = 1;
                    break;
                }
            }
        }
        if(exist_flag == 0) {
            exist_sites[p.getId()].push_back(p.getStart());
            primer_vector_new.push_back(p);
        }
    }
    return primer_vector_new;
}

void generatePrimerPairResult_old(vector<Primer> primer_vector_left, vector<Primer> primer_vector_right, int min_product_size, int max_product_size, int diff_setting, string end_flag, string output, string output_download) {
    ofstream fout(output);
    ofstream fout_download(output_download);
    fout_download << "index" << "\t" << "seq id" << "\t" << "type" << "\t" << "primer start" << "\t" << "primer length" << "\t" << "porduct size" << "\t" << "Tm" << "\t" << "GC content" << "\t" << "any" << "\t" << "3\'" << "\t" << "diff count in primer" << "\t" << "3\' end diff" << "\t" << "primer seq" << "\n";

    vector<Primer_Pair> ppv;
    for(int i = 0; i < primer_vector_left.size(); i++) {
        Primer p_l = primer_vector_left[i];

        if(p_l.getDiff_site_number() < diff_setting) {
            continue;
        }
        if(p_l.getDeletion_site_number() > 10) {
            continue;
        }
        /*
        if(end_flag.compare("Yes") == 0) {
            if(p_l.getDiff_3_end_flag().compare(end_flag) != 0) {
                continue;
            }
        }
        */

        for(int j = 0; j < primer_vector_right.size(); j++) {
            Primer p_r = primer_vector_right[j];

            if(p_r.getDiff_site_number() < diff_setting) {
                continue;
            }
            if(p_r.getDeletion_site_number() > 10) {
                continue;
            }
            /*
            if(end_flag.compare("Yes") == 0) {
                if(p_r.getDiff_3_end_flag().compare(end_flag) != 0) {
                    continue;
                }
            }
            */

            if(p_l.getId().compare(p_r.getId()) == 0) {
                /*
                if(p_l.getId().find("1DL") > p_l.getId().length()) {
                    continue;
                }
                */
                string diff_l_str = generateDifferentSiteStringInPrimer(p_l);
                string diff_r_str = generateDifferentSiteStringInPrimer(p_r);
                int product_size = p_r.getStart() - p_l.getStart() + 1;
                if(product_size >= min_product_size && product_size <= max_product_size) {
                    //primer_pair_count++;
                    //cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                    //fout << primer_pair_count << "\t" + p_l.getId() + "\t" + "LEFT_PRIMER" + "\t" << p_l.getStart() << "\t" << p_l.getLength() << "\t" << product_size << "\t" + p_l.getTm() + "\t" + p_l.getGc() + "\t" + p_l.getAny() + "\t" + p_l.getThree() + "\t" << p_l.getDiff_site_number() << "\t" + p_l.getDiff_3_end_flag() + "\t" + p_l.getSeq() << "\t" << p_l.getStart_in_alignment() << "\t" << p_l.getEnd_in_alignment() << "\t" << diff_l_str << "\n";
                    //fout << primer_pair_count << "\t" + p_r.getId() + "\t" + "RIGHT_PRIMER" + "\t" << p_r.getStart() << "\t" << p_r.getLength() << "\t" << product_size << "\t" + p_r.getTm() + "\t" + p_r.getGc() + "\t" + p_r.getAny() + "\t" + p_r.getThree() + "\t" << p_r.getDiff_site_number() << "\t" + p_r.getDiff_3_end_flag() + "\t" + p_r.getSeq() << "\t" << p_r.getStart_in_alignment() << "\t" << p_r.getEnd_in_alignment() << "\t" << diff_r_str << "\n";

                    //fout_download << primer_pair_count << "\t" + p_l.getId() + "\t" + "LEFT_PRIMER" + "\t" << p_l.getStart() << "\t" << p_l.getLength() << "\t" << product_size << "\t" + p_l.getTm() + "\t" + p_l.getGc() + "\t" + p_l.getAny() + "\t" + p_l.getThree() + "\t" << p_l.getDiff_site_number() << "\t" + p_l.getDiff_3_end_flag() + "\t" + p_l.getSeq() << "\n";
                    //fout_download << primer_pair_count << "\t" + p_r.getId() + "\t" + "RIGHT_PRIMER" + "\t" << p_r.getStart() << "\t" << p_r.getLength() << "\t" << product_size << "\t" + p_r.getTm() + "\t" + p_r.getGc() + "\t" + p_r.getAny() + "\t" + p_r.getThree() + "\t" << p_r.getDiff_site_number() << "\t" + p_r.getDiff_3_end_flag() + "\t" + p_r.getSeq() << "\n";
                    if(end_flag.compare("Yes") == 0) {
                        if(p_l.getDiff_3_end_flag().compare(end_flag) == 0 || p_r.getDiff_3_end_flag().compare(end_flag) == 0) {
                            Primer_Pair pp;
                            pp.setLeft_primer(p_l);
                            pp.setRight_primer(p_r);
                            pp.setScore(p_l.getScore() + p_r.getScore());
                            ppv.push_back(pp);
                        }
                    } else {
                        Primer_Pair pp;
                        pp.setLeft_primer(p_l);
                        pp.setRight_primer(p_r);
                        pp.setScore(p_l.getScore() + p_r.getScore());
                        ppv.push_back(pp);
                    }
                }
            }
        }
    }

    //sort(primer_vector.begin(), primer_vector.end(), cmpScore);
    sort(ppv.begin(),ppv.end(), cmpPairScore);

    int primer_pair_count = 0;
    for(int i = 0; i < ppv.size(); i++) {
        Primer p_l = ppv[i].getLeft_primer();
        Primer p_r = ppv[i].getRight_primer();

        string diff_l_str = generateDifferentSiteStringInPrimer(p_l);
        string diff_r_str = generateDifferentSiteStringInPrimer(p_r);
        int product_size = p_r.getStart() - p_l.getStart() + 1;

        primer_pair_count++;

        fout << primer_pair_count << "\t" + p_l.getId() + "\t" + "LEFT_PRIMER" + "\t" << p_l.getStart() << "\t" << p_l.getLength() << "\t" << product_size << "\t" + p_l.getTm() + "\t" + p_l.getGc() + "\t" + p_l.getAny() + "\t" + p_l.getThree() + "\t" << p_l.getDiff_site_number() << "\t" + p_l.getDiff_3_end_flag() + "\t" + p_l.getSeq() << "\t" << p_l.getStart_in_alignment() << "\t" << p_l.getEnd_in_alignment() << "\t" << diff_l_str << "\n";
        fout << primer_pair_count << "\t" + p_r.getId() + "\t" + "RIGHT_PRIMER" + "\t" << p_r.getStart() << "\t" << p_r.getLength() << "\t" << product_size << "\t" + p_r.getTm() + "\t" + p_r.getGc() + "\t" + p_r.getAny() + "\t" + p_r.getThree() + "\t" << p_r.getDiff_site_number() << "\t" + p_r.getDiff_3_end_flag() + "\t" + p_r.getSeq() << "\t" << p_r.getStart_in_alignment() << "\t" << p_r.getEnd_in_alignment() << "\t" << diff_r_str << "\n";

        fout_download << primer_pair_count << "\t" + p_l.getId() + "\t" + "LEFT_PRIMER" + "\t" << p_l.getStart() << "\t" << p_l.getLength() << "\t" << product_size << "\t" + p_l.getTm() + "\t" + p_l.getGc() + "\t" + p_l.getAny() + "\t" + p_l.getThree() + "\t" << p_l.getDiff_site_number() << "\t" + p_l.getDiff_3_end_flag() + "\t" + p_l.getSeq() << "\n";
        fout_download << primer_pair_count << "\t" + p_r.getId() + "\t" + "RIGHT_PRIMER" + "\t" << p_r.getStart() << "\t" << p_r.getLength() << "\t" << product_size << "\t" + p_r.getTm() + "\t" + p_r.getGc() + "\t" + p_r.getAny() + "\t" + p_r.getThree() + "\t" << p_r.getDiff_site_number() << "\t" + p_r.getDiff_3_end_flag() + "\t" + p_r.getSeq() << "\n";
    }

    fout.close();
    fout_download.close();
}

void generatePrimerPairResult(vector<Primer> primer_vector_left, vector<Primer> primer_vector_right, int min_product_size, int max_product_size, string filter_type, int total_diff, int terminal_diff, int terminal_len, string last_base_diff_flag) {
    int t = 0;
    if(filter_type.compare("AAA") == 0) {
        t = 1;
    }
    if(filter_type.compare("AAB") == 0) {
        t = 2;
    }
    if(filter_type.compare("ABA") == 0) {
        t = 3;
    }
    if(filter_type.compare("ABB") == 0) {
        t = 4;
    }
    if(filter_type.compare("BAA") == 0) {
        t = 5;
    }
    if(filter_type.compare("BAB") == 0) {
        t = 6;
    }
    if(filter_type.compare("BBA") == 0) {
        t = 7;
    }
    if(filter_type.compare("BBB") == 0) {
        t = 8;
    }

    int primer_pair_count = 0;
    for(int i = 0; i < primer_vector_left.size(); i++) {
        Primer p_l = primer_vector_left[i];
        for(int j = 0; j < primer_vector_right.size(); j++) {
            Primer p_r = primer_vector_right[j];
            if(p_l.getId().compare(p_r.getId()) == 0) {
                if(p_l.getId().find("1BL") > p_l.getId().length()) {
                    continue;
                }
                if(p_r.getStart() - p_l.getStart() + 1 >= min_product_size && p_r.getStart() - p_l.getStart() + 1 <= max_product_size) {
                    int tmp_l = calculateDiffSiteIn3Terminal(p_l, terminal_len, "LEFT_PRIMER");
                    int tmp_r = calculateDiffSiteIn3Terminal(p_r, terminal_len, "RIGHT_PRIMER");
                    switch(t) {
                        case 1:
                            if((p_l.getDiff_site_number() >= total_diff && p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff  && tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 && p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 2:
                            if((p_l.getDiff_site_number() >= total_diff && p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff  && tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 || p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 3:
                            if((p_l.getDiff_site_number() >= total_diff && p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff || tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 && p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 4:
                            if((p_l.getDiff_site_number() >= total_diff && p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff || tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 || p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 5:
                            if((p_l.getDiff_site_number() >= total_diff || p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff  && tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 && p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 6:
                            if((p_l.getDiff_site_number() >= total_diff || p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff && tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 || p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 7:
                            if((p_l.getDiff_site_number() >= total_diff || p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff || tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 && p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                        case 8:
                            if((p_l.getDiff_site_number() >= total_diff || p_r.getDiff_site_number() >= total_diff) && (tmp_l >= terminal_diff || tmp_r >= terminal_diff) && (p_l.getDiff_3_end_flag().compare(last_base_diff_flag) == 0 || p_r.getDiff_3_end_flag().compare(last_base_diff_flag) == 0)) {
                                primer_pair_count++;
                                cout << p_l.getId() << " " << p_l.getStart() << " " << p_l.getSeq() << " " << p_l.getDiff_site_number() << " " << p_l.getDiff_3_end_flag() << " " << p_l.getScore() << " " << p_r.getId() << " " << p_r.getStart() << " " << p_r.getSeq() << " " << p_r.getDiff_site_number() << " " << p_r.getDiff_3_end_flag() << " " << p_r.getScore() << "\n";
                            }
                            break;
                    }
                }
            }
        }
    }
}

int getAlignmentSite(string alignment_seq, int real_site) {
    int gap_base_number = 0;
    int no_gap_base_number= 0;
    int alignmnet_site = 0;
    for(int i = 0; i < alignment_seq.length(); i++) {
        if(alignment_seq.at(i) == '-') {
            gap_base_number++;
        } else {
            no_gap_base_number++;
        }
        if(no_gap_base_number == real_site) {
            alignmnet_site = i + 1;
            break;
        }
    }
    return alignmnet_site;
}

string generateDifferentSiteStringInPrimer(Primer p) {
    vector<DifferenceSite> tmp_v = p.getDiff_vector();
    string str = "";
    for(int i = 0; i < tmp_v.size(); i++) {
        str += boost::lexical_cast<string>(tmp_v[i].getIn_alignment_site()) + ";";
    }
    return str.substr(0, str.length() - 1);
}

string getReverseComplement(string seq) {
    string new_seq;
    for(int i = seq.length() - 1; i >= 0; i--) {
        char c = seq.at(i);
        if(c == 'A' || c == 'a') {
            new_seq += "T";
        } else if (c == 'T' || c == 't') {
            new_seq += "A";
        } else if (c == 'C' || c == 'c') {
            new_seq += "G";
        } else if (c == 'G' || c == 'g') {
            new_seq += "C";
        } else {
            new_seq += c;
        }
    }
    return new_seq;
}

map<string, vector<DifferenceSite>> identifyDifferentSite(map<string, string> alignment_seq_map) {
    map<string, string>::iterator iter = alignment_seq_map.begin();
    int length = 0;
    vector<string> ids;
    for(; iter != alignment_seq_map.end(); ++iter) {
        length = iter->second.length();
        ids.push_back(iter->first);
    }

    map<string, vector<DifferenceSite>> differentsite_map;
    for(int i = 0; i < length; i++) {
        map<string, string>::iterator it = alignment_seq_map.begin();
        string consesus;
        for(; it != alignment_seq_map.end(); ++it) {
            char s = it->second.at(i);
            consesus += s;
        }

        for(int j = 0; j < consesus.length(); j++) {
            int char_count = 0;
            char a = consesus.at(j);
            for(int k = 0; k < consesus.length(); k++) {
                if(a == consesus.at(k)) {
                    char_count++;
                }
            }
            if(char_count == 1) {
                //cout << consesus << " " << ids[j] << " " << i + 1 << "\n";
                if(differentsite_map.find(ids[j]) != differentsite_map.end()) {
                    vector<DifferenceSite> tmp_v = differentsite_map[ids[j]];
                    DifferenceSite ds;
                    ds.setId(ids[j]);
                    ds.setIn_alignment_site(i + 1);
                    string tmp_seq = alignment_seq_map[ids[j]];
                    ds.setIn_seq_site(calculateSeqSitebyAlignmentSite(tmp_seq, i+ 1));
                    ds.setChar_in_different(a);
                    tmp_v.push_back(ds);
                    differentsite_map[ids[j]] = tmp_v;
                } else {
                    vector<DifferenceSite> tmp_v;
                    DifferenceSite ds;
                    ds.setId(ids[j]);
                    ds.setIn_alignment_site(i + 1);
                    string tmp_seq = alignment_seq_map[ids[j]];
                    ds.setIn_seq_site(calculateSeqSitebyAlignmentSite(tmp_seq, i+ 1));
                    ds.setChar_in_different(a);
                    tmp_v.push_back(ds);
                    differentsite_map[ids[j]] = tmp_v;
                }
            }
        }
    }
    //cout << differentsite_map["chr1AL"].size() << "\n";
    return differentsite_map;
}

int calculateSeqSitebyAlignmentSite(string seq, int alignment_site) {
    int count = 0;
    for(int i = 0; i < seq.length(); i++) {
        if(i + 1 > alignment_site) {
            break;
        }
        char a = seq.at(i);
        if(a == '-') {
            count++;
        }
    }
    return alignment_site - count;
}

int calculateDiffSiteIn3Terminal(Primer p, int len, string type) {
    int count_diff_in_3_terminal = 0;
    vector<DifferenceSite> diff_vector = p.getDiff_vector();
    //cout << diff_vector.size() << "\n";
    for(int i = 0; i < diff_vector.size(); i++) {
        DifferenceSite ds = diff_vector[i];
        if(type.compare("LEFT_PRIMER") == 0) {
            if(ds.getIn_seq_site() >= p.getStart() + p.getLength() - len && ds.getIn_seq_site() < p.getStart() + p.getLength()) {
                count_diff_in_3_terminal++;
            }
        } else {
            if(ds.getIn_seq_site() > p.getStart() - p.getLength() && ds.getIn_seq_site() <= p.getStart() - p.getLength() + len) {
                count_diff_in_3_terminal++;
            }
        }
    }
    return count_diff_in_3_terminal;
}

void calculatePrimerDimer(string primer3_output) {
    ifstream fin(primer3_output);
    string line;
    int flag = 0;
    string seq_id;
    ofstream fout_l;
    ofstream fout_r;
    while(getline(fin, line)) {
        if(line.find("PRIMER PICKING RESULTS FOR ") < line.length()) {
            if(fout_l) {
                fout_l.close();
            }
            if(fout_r) {
                fout_r.close();
            }
            flag = 1;
            seq_id = replace_all(line, "PRIMER PICKING RESULTS FOR ", "");
            fout_l.open(primer3_output.substr(0, primer3_output.find_last_of("/")) + "/" + seq_id + "_l");
            fout_r.open(primer3_output.substr(0, primer3_output.find_last_of("/")) + "/" + seq_id + "_r");
        }
        if(flag == 1 && line.find(" LEFT_PRIMER") < line.length()) {
            boost::algorithm::trim(line);
            boost::xpressive::sregex rgx = sregex::compile("\\s+");
            boost::xpressive::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
            boost::xpressive::sregex_token_iterator end;

            int index = 0;
            while(iter != end) {
                //cout << *iter << "\n";
                if(index == 0) {
                    string index_id = *iter;
                    fout_l << ">" + seq_id + "__" + index_id << "\n";
                } else if(index == 8) {
                    string seq = *iter;
                    fout_l << seq << "\n";
                    break;
                }

                index++;
                iter++;
            }
        }
        if(flag == 1 && line.find(" RIGHT_PRIMER") < line.length()) {
            boost::algorithm::trim(line);
            boost::xpressive::sregex rgx = sregex::compile("\\s+");
            boost::xpressive::sregex_token_iterator iter(line.begin(), line.end(), rgx, -1);
            boost::xpressive::sregex_token_iterator end;

            int index = 0;
            while(iter != end) {
                //cout << *iter << "\n";
                if(index == 0) {
                    string index_id = *iter;
                    fout_r << ">" + seq_id + "__" + index_id << "\n";
                } else if(index == 8) {
                    string seq = *iter;
                    fout_r << seq << "\n";
                    break;
                }

                index++;
                iter++;
            }
        }
    }
    if(fout_l) {
        fout_l.close();
    }
    if(fout_r) {
        fout_r.close();
    }
    fin.close();
}

string& replace_all(string&   str,const   string&   old_value,const   string&   new_value)
{
    while(true)   {
        string::size_type   pos(0);
        if(   (pos=str.find(old_value))!=string::npos   )
            str.replace(pos,old_value.length(),new_value);
        else   break;
    }
    return   str;
}