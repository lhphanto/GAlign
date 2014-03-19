// class of global alignment algorithm of two vectors of strings
#ifndef _ALIGN_H_
#define _ALIGN_H_

#include<string>
#include<vector>
#include<iostream>
#include<math.h>
#include<map>

using std::map;

using namespace std;


class Alignment{
	private:
		static const MaxSeq = 4; //for safety reason we now align at most 4 sequences
		//store of the two original sequences
		vector<string> seq1;
		vector<string> seq2;
		//store of the aligned sequences
		vector<string> seq1a;
		vector<string> seq2a;
		string name1;
		string name2;
		//penalty scores for gap,mismatch and similar
		float gap;
		float mismatch;
		float similar;
		float identical;
		float match;
		bool FloatEqual(float A, float B);
	public:
		Alignment(vector<string>& a,vector<string>& b);
		Alignment(vector<string>& a,string namea,vector<string>& b,string nameb);
		Alignment(vector<string>& a,vector<string>& b, float,float,float,float,float);

		~Alignment();		
		//Main function to do the alignment
		//Store result in seq1a,seq2a
		//Return the normalized score of the alignment
		float align();

		//print out the two aligned sequences
		void outputalign();

		//store the seq1a into 'store' if index == 1, store the seq2a if index == 2  
		void outputalign(vector<string> & store, int index);

		//Input two strings
		//Output fraction of the characters that are the same
		float Similarity(string a,string b);
};

#endif 
