#include<cstdlib>
#include<fstream>
#include "align.h"
#include "malign.h"

using namespace std;

void Tokenize(const string& , vector<string>&, const string&);

int main(int argc, char* argv[]){
	//cout << "argc = " << argc << endl; 
	ifstream inFile;
	string line;
	vector<string> codes1;
	vector<string> codes2;
	vector<string> proteins;
	vector<string> OrderPro;
	float score;
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
			proteins.push_back(line);
		}
	}
	OrderPro.push_back(proteins[0]);
	proteins.erase (proteins.begin());
	int count=0;
	while(proteins.size() > 0){
		float minScore=100.0;
		int minIndex;
		for(unsigned int i=0; i < proteins.size();i++){
			codes1.clear();
			codes2.clear();
			Tokenize(OrderPro[OrderPro.size()-1],codes1,"\t");
			Tokenize(proteins[i],codes2,"\t");
			if(codes1.size() < 3 || codes2.size() < 3){
				cout<<"Error:Either the format is not tab-delimited or the code sequence is too short"<<endl;
				return 1;
			}	
			Alignment Temp(codes1,codes2);
			score=Temp.align();
			if(minScore == 100.0 || minScore >score){
				minScore=score;
				minIndex=i;
			}
		}
		OrderPro.push_back(proteins[minIndex]);
		proteins.erase (proteins.begin()+minIndex);
		count++;
	//	cout<<"Index:"<<count<<",Score:"<<score<<endl;
			//~ Temp.outputalign();
			
		//~ for(int j=i+1;j <proteins.size();j++){
			
		//~ }
	}
	//cout<<"Finish Pairwise?"<<endl;
	MAlignment X;
	for(unsigned int i=0; i < OrderPro.size();i++){
		codes1.clear();
		Tokenize(OrderPro[i],codes1,"\t");
		X.AddSeq(codes1);
	}
	float myscore=X.align();
	//myscore = myscore/(float)OrderPro.size();
	cout<<myscore<<endl;
	X.outputalign();	
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


