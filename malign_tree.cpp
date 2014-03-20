#include "malign_tree.h"

MalignmentT::MalignmentT(vector< vector<string> > input_data){
	for(unsigned int i=0; i < input_data.size();i++){
		string name = input_data[i][0];
		input_data[i].erase(input_data[i].begin());
		this->seq_data.push_back(input_data[i]);
		this->ind2name.insert(std::pair<int,string>(i,name));
	}
	this->Num_Seq = (int)input_data.size();
	
	this->dis_mtx = new double[this->Num_Seq*this->Num_Seq];
	for(int i=0; i < this->Num_Seq;i++)
		for(int j=0; j < this->Num_Seq;j++)
			this->dis_mtx[i*this->Num_Seq+j]=100.0;

	//use the same setting as in malign.h
	this->gap = 1.0;
	this->mismatch = 1.0;
	this->similar = -0.5;
	this->verysimilar = -1.0;
	this->identical = -2.0; 
	return;
}		

MalignmentT::~MalignmentT(){
	delete[] this->dis_mtx;
	this->dis_mtx = NULL;
}

int MalignmentT::align(void){
	cal_dismtx();

	/*for testing cal_upgma tree	
	double debug[49] = {0.0,19.0,27.0,8.0,33.0,18.0,13.0,
				19.0,0,31.0,18.0,36.0,1.0,13.0,
				27.0,31.0,0,26.0,41.0,32.0,29.0,
				8.0,18.0,26.0,0,31.0,17.0,14.0,
				33.0,36.0,41.0,31.0,0,35.0,28.0,
				18.0,1.0,32.0,17.0,35.0,0,12.0,
				13.0,13.0,29.0,14.0,28.0,12.0,0};
	
	cal_upgma(debug,7);
	*/
	cal_upgma(this->dis_mtx,this->Num_Seq,true);
	//print_tree();
	//print_mtx();
	//print_alignment();
	return 0;
}


//aligning the alignment profiles in b and c and store
//the new profile into a
//the table is like :
//------------b(x)-------
//|
//|
//|
//a 
//|
//|
//trace = 1 (going right), 3(going down),2(going diagonally)
float MalignmentT::align_prof(vector< vector<string> >& a, vector< vector<string> >& b, vector< vector<string> >& c){
	
	float FinalScore;

	
	/*******Calculate the TABLE********/		
	int ydim=(int)a[0].size()+1;
	int xdim=(int)b[0].size()+1;
	float score[xdim*ydim];
	int trace[xdim*ydim]; // to record traces for finding the aligned sequences 1 for gap in y,2 for match,3 for gap in x
	for(int l=0; l < ydim;l++){
		score[l*xdim]=(float)l*gap;
		trace[l*xdim]=3;
	}
	for(int l=0; l < xdim;l++){
		score[l]=(float)l*gap;
		trace[l]=1;
	}

	for(int j=1; j <xdim;j++){
		for(int k=1; k <ydim;k++){
			float  temp=0.0;
			score[k*xdim+j]=score[(k-1)*xdim+j]+gap;//initialize min with matching y dim to gap
			trace[k*xdim+j]=3;
			
			for(unsigned int m=0;m < b.size();m++){
				for(unsigned int n=0;n < a.size();n++){
					//cout<<"DEBUG:"<<b[0].size()<<","<<k<<":"<<a[0].size()<<"."<<endl;
					temp= temp + Score(b[m][j-1],a[n][k-1]);
				}
			}
			temp=temp/float(b.size()*a.size());
			temp=temp+score[(k-1)*xdim+j-1];
			
			if(temp < score[k*xdim+j]){
				
				score[k*xdim+j]=temp;
				trace[k*xdim+j]=2;
			}

			temp=score[k*xdim+j-1]+gap; //adding gap in x direction
			if(temp < score[k*xdim+j]){
				//cout<<temp << " is smaller than "<< score[k*xdim+j] <<endl;
				score[k*xdim+j]=temp;
				trace[k*xdim+j]=1;
			}
		}
	}
		//cout<<"Finish Scoring?"<<endl;
		
		/**************  Output score matrix*************
		for(int k=ydim-1; k >=0;k--){
			if( k >= 1){
				cout<<profile[0][k-1];
			}else{
				cout<<"-";
			}
			for(int j=0; j <xdim;j++){
				cout<<"\t"<<score[k*xdim+j];
			}
			cout<<endl;
		}
		cout<<"\t"<<"-";
		for(unsigned int j=0; j < seqs[i].size();j++){
			cout<<"\t"<<seqs[i][j];
		}
		cout<<endl;
		cout<<"Final Score :"<<score[xdim*ydim-1]<<endl;
		*************************************************/
		
		/**************  Out put matrix
		for(int l=ydim-1; l >=0;l--){			
			for(int m=0; m <xdim;m++){
				cout<<"\t"<<trace[l*xdim+m];
			}
			cout<<endl;
		}
		*/
		
	/**Trace back find the alignment*/
	int xInd=xdim-1;
	int yInd=ydim-1;

	//renew c using a
	c = a;	
	for(unsigned int i=0; i < b.size();i++)
		c.push_back(b[i]);
	while(xInd > 0 || yInd > 0 ){
		//cout <<"Checking "<<xInd<<","<<yInd<<endl;
		//gap in x(b) axis, 
		if(yInd > 0 && trace[yInd*xdim+xInd] == 3){
			for(unsigned int m=a.size();m < c.size();m++)
				c[m].insert(c[m].begin()+xInd,"-");	
			yInd=yInd-1;
			continue;
		}
			
		if(xInd > 0 && trace[yInd*xdim+xInd] == 1){
			//cout<<xInd<<endl;
			for(unsigned int m=0;m < a.size();m++)
				c[m].insert(c[m].begin()+yInd,"-");					
			xInd=xInd-1;
			continue;
		}
		if(trace[yInd*xdim+xInd] == 2){				
			xInd=xInd-1;
			yInd=yInd-1;
		}
	}
		
	FinalScore=score[xdim*ydim-1];
	
	return (FinalScore)/(float)(c[0].size())+2.0;

}

void MalignmentT::cal_dismtx(void){
	double score;
	
	for(int i=0; i < this->Num_Seq ;i++){
		for(int j=i; j < this->Num_Seq;j++){
			//Alignment Temp(this->seq_data[i],"a",this->seq_data[j],"b");
			Alignment Temp;
			Temp.AddSeq(this->seq_data[i],"a");
			Temp.AddSeq(this->seq_data[j],"b");
			//cout<<this->ind2name[i]<<"=>"<<this->seq_data[i].size()<<endl;
			score=Temp.Align();
			//cout<<"Score is "<<score<<endl;	
			this->dis_mtx[i*this->Num_Seq+j]=score;
			this->dis_mtx[j*this->Num_Seq+i]=score;
			
		}
	}


	for(int i=0; i < this->Num_Seq ;i++){
		double shift = this->dis_mtx[i*this->Num_Seq+i];
		for(int j=0; j < this->Num_Seq;j++)
			this->dis_mtx[i*this->Num_Seq+j]=this->dis_mtx[i*this->Num_Seq+j]-shift;
	}

	
	return;	
}

//calculate score for match of a to b
//a b must be of equal length
float  MalignmentT::Score(string a,string b){
	//cout<<"Comparing "<<a<<","<<b<<endl;
	if(a.compare("-") == 0 || b.compare("-") == 0)
		return this->gap;
	if(a.size() != b.size()){
		cerr<<"comparing string elements of different length:"<<a<<","<<b<<endl;
		exit(1);
	}

	float sim=0.0;
	string tempa,tempb;
	for(unsigned int i=0;i < a.size();i++){
		tempa=a.substr(i,1);
		tempb=b.substr(i,1);
		if(tempa.compare(tempb)== 0 && tempa.compare("x") != 0 ){
			sim=sim+1.0;
		}else{
			if(a.length() == 4 && i == 1){
				sim=sim+0.5;
			}
		}
	}

	if( sim == ((float)a.size()-0.5) )
		return this->verysimilar;

	if( sim == (float)a.size() )
		return this->identical;
	
	if( sim > (float)a.size()*0.66 )
		return this->similar;
		
	return this->mismatch;	
}

void MalignmentT::find_min_pair(double& min_dist,int& a,int& b){
	min_dist = this->guide_tree[0]->dists[1];
	a =0;
	b =1;
	for(unsigned int i=0; i < this->guide_tree.size();i++){
		for(unsigned int j=i+1; j < this->guide_tree[i]->dists.size();j++){
			if( this->guide_tree[i]->dists[j] < min_dist){
				min_dist = this->guide_tree[i]->dists[j];
				a = i;
				b = j;
			}
		}
	}
	return;
}		

void MalignmentT::update_dist(int a,int b){
	for(unsigned int i=0; i < this->guide_tree.size();i++){
		if((int)i == a || (int)i == b)
			continue;
		double sizea = (double)this->guide_tree[a]->size;
		double sizeb = (double)this->guide_tree[b]->size;
		double new_dist = (sizea*this->guide_tree[i]->dists[a]+sizeb*this->guide_tree[i]->dists[b])/(sizea+sizeb); 			
		this->guide_tree[i]->dists.erase(this->guide_tree[i]->dists.begin()+a);	
		this->guide_tree[i]->dists.erase(this->guide_tree[i]->dists.begin()+b-1);	
		this->guide_tree[i]->dists.push_back(new_dist);
	}
	return;
}

void MalignmentT::cal_upgma(double* dist_mtx,int dim,bool align){
	//Initialize the tree nodes
	for(int i=0; i < dim;i++){
		TreeNode* temp = new TreeNode();	
		temp->seq_ind = i;
		temp->align_profile.push_back(this->seq_data[i]);
		temp->align_seq_inds.push_back(i);
		for(int j=0;j< dim;j++)
			temp->dists.push_back(dist_mtx[i*dim+j]);
		this->guide_tree.push_back(temp);
	}

	//converging and construct the tree
	while(this->guide_tree.size() > 1){
		double min_dist = -100.0;
		int i,j;
		find_min_pair(min_dist,i,j);
		update_dist(i,j);
		
		TreeNode* temp = new TreeNode();
		temp->size = this->guide_tree[i]->size+this->guide_tree[j]->size;
		temp->depth = min_dist/2.0;
		temp->lchild = this->guide_tree[i];
		temp->rchild = this->guide_tree[j];
		temp->ldist = temp->depth - temp->lchild->depth;
		temp->rdist = temp->depth - temp->rchild->depth; 
		this->guide_tree.erase(this->guide_tree.begin()+i);
		this->guide_tree.erase(this->guide_tree.begin()+j-1);
		for(unsigned int k=0; k < this->guide_tree.size();k++){
			int array_size = (int)this->guide_tree[k]->dists.size();
			//cout<<array_size<<","<<this->guide_tree.size()<<endl;
			temp->dists.push_back(this->guide_tree[k]->dists[array_size-1]);
		}
		temp->dists.push_back(0.0);
		if(align == true){
			temp->align_score=align_prof(temp->lchild->align_profile,temp->rchild->align_profile,temp->align_profile);
			temp->align_seq_inds = temp->lchild->align_seq_inds;
			for(unsigned int i=0; i < temp->rchild->align_seq_inds.size();i++)
				temp->align_seq_inds.push_back(temp->rchild->align_seq_inds[i]);
		}
		this->guide_tree.push_back(temp);
		//cout<<"Now the tree is "<<endl;
		//print_tree();
	}	
		
}
void MalignmentT::print_data(void){
	for(map<int,string >::iterator it=this->ind2name.begin(); it != this->ind2name.end(); ++it){
		cout << it->second << "=>";
		for(unsigned int i=0; i < this->seq_data[it->first].size();i++)
			cout<<","<<this->seq_data[it->first][i];
		cout<<endl;
	}
}

void MalignmentT::print_mtx(void){
	for(int i=0; i < this->Num_Seq ;i++){
		for(int j=0; j < this->Num_Seq;j++)
			cout<<this->dis_mtx[i*this->Num_Seq+j]<<" ";
		cout<<endl;
	}
	return;
}

void MalignmentT::print_treenode(TreeNode* node){
	if(node == NULL)
		return;
		
	if(node->seq_ind >= 0){
		cout<<this->ind2name[node->seq_ind];
	}else{
		cout<<"(";
		print_treenode(node->lchild);
		cout<<":"<<node->ldist;
		cout<<",";
		print_treenode(node->rchild);
		cout<<":"<<node->rdist<<")";
	}
	return;
}


void MalignmentT::print_tree(void){
	for(unsigned int i=0; i < this->guide_tree.size();i++){
		print_treenode(this->guide_tree[i]);
		cout<<endl;
	}
}

void MalignmentT::print_alignment(void){
	for(unsigned int i=0; i < this->guide_tree.size();i++){
		cout<<this->guide_tree[i]->align_score<<endl;
		for(unsigned int j=0; j < this->guide_tree[i]->align_profile.size();j++){
			cout<<this->ind2name[this->guide_tree[i]->align_seq_inds[j]];
			for(unsigned int k=0; k < this->guide_tree[i]->align_profile[j].size();k++)
				cout<<"\t"<<this->guide_tree[i]->align_profile[j][k];
			cout<<endl;
		}
	}
	return;
}
	
