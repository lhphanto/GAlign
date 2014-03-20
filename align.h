// class of global alignment algorithm of two vectors of strings
#ifndef _ALIGN_H_
#define _ALIGN_H_

#include<string>
#include<vector>
#include<iostream>
#include<math.h>
#include<bitset>

using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;

struct MtxE{
	MtxE():val(-1),prev(-1){}
	float val;
	int prev;
};

class Alignment{
	private:
		static const int MaxSeq = 9; //for safety reason we now align at most 4 sequences
		//store of the two original sequences
		//also we make sure the first seq is always the shortest one
		vector< vector<string> > inseqs;
		//store of the aligned sequences
		vector< vector<string> > otseqs;
		vector<string> names;
		float FinalScore; 
		//penalty scores for gap,mismatch and similar
		float Gap;
		float Mismatch;
		float Similar;
		float Identical;
		float Match;
		bool FloatEqual(float A, float B);
		int cor_to_ind(int*);//translate cors to index in Mtx
		void reset_cor(int*);//reset cors array to 0
		float pos_score(MtxE*,int*,bool*);
		void RecCal(int lev,int* Cors,MtxE* mtx,vector<int> nzero);
		void TraceAln(MtxE*);
		void outputScoreMtxHelper(int,int*,MtxE*);
		void outputScoreMtx(MtxE* mtx);
	public:
		Alignment();
		//initialize the alignment by specify your own scores
		Alignment(float,float,float,float,float);

		~Alignment();
		bool AddSeq(vector<string> const &);	
		bool AddSeq(vector<string> const &,string const &);	
		bool AddSeqs(vector< vector<string> > const &);	
		bool AddSeqs(vector< vector<string> > const &,vector<string> const &);	


		//Main function to do the alignment
		//Store result in seq1a,seq2a
		//Return the normalized score of the alignment
		float Align();

		//print out the two aligned sequences
		void PrintAlignment();

		//store the seq1a into 'store' if index == 1, store the seq2a if index == 2  
		void PrintAlignment(vector<string> & store, int index);

		//Input two strings
		//Output fraction of the characters that are the same
		float Similarity(string a,string b);
};

#endif 
