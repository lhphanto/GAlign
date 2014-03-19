#include<cstdlib>
#include<fstream>
#include "align.h"
#include "malign_tree.h"

using namespace std;

void Tokenize(const string& , vector<string>&, const string&);

int main(int argc, char* argv[]){

	//cout << "argc = " << argc << endl; 
	ifstream inFile;
	string line;
	vector<	vector<string> > data_mtx;
	vector<string> codes;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open file datafile.txt";
		exit(1);   // call system to stop
	}
	
	
	
	//pairwise alignment
	while(!inFile.eof()){
		getline(inFile,line);
		if(line.size() > 0){
			codes.clear();
			Tokenize(line,codes,"\t");
			data_mtx.push_back(codes);	
		}
	}
	
	MalignmentT aln(data_mtx);
	//aln.print_data();
	aln.align();
	aln.print_alignment();
	return 0;
}


void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}


