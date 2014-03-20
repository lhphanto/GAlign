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
	vector<string> codes;
	vector< vector<string> >data_mtx;
	//int GraphSize=0;
	inFile.open(argv[1]);
	if (!inFile) {
		cerr << "Unable to open file datafile.txt";
		exit(1);   // call system to stop
	}
		
	//~ getline(inFile,line);	
	//~ Tokenize(line,codes1,"\t");	
	//~ getline(inFile,line);	
	//~ Tokenize(line,codes2,"\t");
	
		//~ //cout << genes[1]<<"."<<genes.size()<<endl;
	//~ if(codes1.size() >0 && codes2.size() >0){
		//~ cout<<codes1[0]<<endl;
		//~ cout<<codes2[0]<<endl;
		//~ codes1.erase (codes1.begin());
		//~ codes2.erase (codes2.begin());
		//~ Alignment A(codes1,codes2,3,2);
		//~ score=A.align();
		//~ cout<<"Score:"<<score<<endl;
	//~ }	
	
	
	//pairwise alignment
	while(!inFile.eof()){
		getline(inFile,line);
		if(line.size() > 0){
			codes.clear();
			Tokenize(line,codes,"\t");
			data_mtx.push_back(codes);
		}
	}

	Alignment X;
	X.AddSeqs(data_mtx);
	X.Align();
	X.PrintAlignment();	
	/*
	MalignmentT X(data_mtx);
	float myscore=X.align();
	//myscore = myscore/(float)OrderPro.size();
	cout<<myscore<<endl;
	X.print_alignment();	
	*/
	//~ for(int i=0; i < OrderPro.size();i++){
		//~ cout<<OrderPro[i]<<endl;
	//~ }
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


