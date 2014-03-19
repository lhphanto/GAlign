#ifndef _MALIGN_TREE_H_
#define _MALIGN_TREE_H_

#include <string>
#include <vector>
#include <cstdlib>
#include <map>
#include "align.h"

using std::string;
using std::map;
using std::vector;
using std::cout;
using std::endl;
using std::cerr;

class MalignmentT{
	private:
		vector< vector<string> > seq_data;
		map<int,string> ind2name;
		double* dis_mtx;
		int Num_Seq;
		struct TreeNode
		{
			TreeNode* lchild;
			TreeNode* rchild;
			double ldist;
			double rdist;
			double depth;//height of this subtree,used for calculating ldist and rdist	
			vector<double> dists; //node distance to other nodes, only needed for tree root nodes
			int seq_ind; //only leaf node have valid seq_ind, which is the index in seq_data
			int size; //size of the tree OR number of leaf nodes OR number of sequences contained
					
			//these two are to store the alignment and the sequence index of the contained sequences
			vector< vector<string> > align_profile;
			vector<int> align_seq_inds;
			float align_score;
			TreeNode():lchild(NULL),rchild(NULL),ldist(-1.0),rdist(-1.0),depth(0.0),seq_ind(-1),size(1){}
		};			
		vector<TreeNode*> guide_tree;
		float gap;
		float mismatch;
		float similar;
		float verysimilar;
		float identical;


		void print_treenode(TreeNode*);
		void find_min_pair(double&,int&,int&);
		void update_dist(int,int);
		float Score(string,string);
		void cal_dismtx(void);
		//do alignment when building tree if boolean tag is set to true
		void cal_upgma(double*,int,bool);
		//aligning the alignment profiles in b and c and store
		//the new profile into a
		float align_prof(vector< vector<string> >& a, vector< vector<string> >& b,vector< vector<string> >& c);
	public:
		MalignmentT(){};
		MalignmentT( vector< vector<string> > );
		~MalignmentT();
		int align(void);
		void print_data(void); //purly for debugging
		void print_mtx(void);
		void print_tree(void);
		void print_alignment(void);
};

#endif	
