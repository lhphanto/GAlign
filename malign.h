#ifndef _MALIGN_H_
#define _MALIGN_H_

#include<string>
#include<vector>
#include<iostream>
#include "align.h"


using namespace std;


class MAlignment{
	private:
		vector<vector<string> > seqs; // This is the original sequences
		vector<vector<string> > profile;//This is the set of aligned sequences
		vector<vector<string> > alignment;//Final alignment
		vector<string> names;
		float gap;
		float mismatch;
		float similar;
		float verysimilar;
		float identical;
		void UpdateProfile();
		bool FloatEqual(float A, float B);
	public:
		//default constructor		
		MAlignment(){
			gap=1.0;
			mismatch=1.0;
			similar=-0.5;
			verysimilar = -1.0;
			identical=-2.0;
		}
		MAlignment(float a,float b,float c);
		void AddSeq(vector<string>& a); // add a single sequence
		void AddSeq(vector<vector<string> > &  a);// add a set of sequences
		float align();//Main program to do the alignment
		//float subalign(int m,int n);
		void outputalign();
		float Similarity(string a,string b);
};

#endif 
