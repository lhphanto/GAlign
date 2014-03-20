#include "align.h"

Alignment::Alignment(){	
	Gap=1.0;
	Mismatch=1.0;
	Identical=-1.0;
	Similar= -0.5;
	Match= -2.0;
	this->MyProSize=1;
}

Alignment::Alignment(float c,float d,float e,float f,float g){
	Gap=c;
	Mismatch=d;
	Identical=e;
	Similar=f;
	Match=g;
	this->MyProSize=1;
}

Alignment::~Alignment(){
}

bool Alignment::AddSeq(vector<string> const & input){
	if((this->MyProSize*((int)input.size()-1)) >= this->MaxProSize){
		cerr<<"MaxSeq limit reached , can't add any more seq"<<endl;
		return 0;
	}
	if(input.size() < 2){
		cerr<<"Seq to add too short, it needs to contain at least a name and a single string "<<endl;
		return 0;
	}
	

	vector<string> newseq = input;
	this->names.push_back(newseq.front());
	newseq.erase(newseq.begin());
	this->inseqs.push_back(newseq);
	this->MyProSize = this->MyProSize*((int)newseq.size());	
	return 1;
}

bool Alignment::AddSeq(vector<string> const & input,string const & iname){
	if((this->MyProSize*((int)input.size())) >= this->MaxProSize){
		cerr<<"MaxSeq limit reached , can't add any more seq"<<endl;
		return 0;
	}
	
	if(input.size() < 1){
		cerr<<"Seq to add too short, it needs to contain at least a single string "<<endl;
		return 0;
	}
	if(iname.size() == 0){
		cerr<<"Input seq's name is an empty string"<<endl;
		return 0;
	}
	
	vector<string> newseq = input;
	this->names.push_back(iname);
	this->inseqs.push_back(newseq);
	this->MyProSize = this->MyProSize*((int)newseq.size());	
	return 1;
}

bool Alignment::AddSeqs(vector< vector<string> > const & input){
	for(unsigned int i=0; i < input.size();i++)
		if(AddSeq(input[i]) == 0)
			return 0;
	return 1;
}

bool Alignment::AddSeqs(vector< vector<string> > const & input,vector<string> const & names){
	if(input.size() != names.size()){
		cerr<<"number of input seqs are not the same as the number of seq names"<<endl;
		return 0;
	}
	for(unsigned int i=0; i < input.size();i++)
		if(AddSeq(input[i],names[i]) == 0)
			return 0;
	return 1;
}



//Helper function to compare two float numbers to 0.01 precision
bool Alignment::FloatEqual(float A, float B){
    if (A == B)
        return true;

    float Error = fabs((A - B));
    if (Error <= 0.01)
        return true;
    return false;

}

//Input two strings
//Output fraction of the characters that are the same
float Alignment::Similarity(string a,string b){
	string tempa,tempb;
	float sim=0.0;
	if(a.length() != b.length() ){
		cerr<<"comparing "<<a<<","<<b<<endl;
		return sim;
	}

	if(a.compare(b) == 0)
		return 1;
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
		return 0.9;

	if( sim == (float)a.size() )
		return 1;
	
		
	return ((float)((int)sim))/(float)a.size();
}

int Alignment::cor_to_ind(int* val){
	int index=0;
	int dim_size = 1;
	int len = (int)this->inseqs.size();
		
	for(int i=len-1;i >=0;i--){
		if(val[i] != 0)
			index += val[i]*dim_size;
		dim_size = dim_size*((int)this->inseqs[i].size()+1);
	}
	return index;
}

void Alignment::reset_cor(int* val){
	for(int i=0;i < (int)this->inseqs.size();i++){
		val[i]=0;
	}	
}
		
float Alignment::Align(void){
	
	if((int)this->inseqs.size() < 2){
		cerr<<"There is only 1 seq instored,can't do alignment"<<endl;
		return 10;
	}
	int TotalSize=1;
	for(unsigned int i=0;i < this->inseqs.size();i++)
		TotalSize = TotalSize*((int)this->inseqs[i].size()+1);

	MtxE * ScoreMtx = new MtxE[TotalSize];
	int* Coors = new int[(int)this->inseqs.size()];
	reset_cor(Coors);
	
	ScoreMtx[0].val=0;
	ScoreMtx[0].prev=-2;
	for(unsigned int i=0; i < this->inseqs.size();i++){
		for(int j=1; j <= (int)this->inseqs[i].size();j++){
			if(i > 0) Coors[i-1]=0;
			Coors[i]=j;
			ScoreMtx[cor_to_ind(Coors)].val=(float)(j*this->Gap);
			ScoreMtx[cor_to_ind(Coors)].prev = 1<<((int)this->inseqs.size()-i-1);//its previous should be(0,0,...,j-1,0,0)	
		}
	}	
	//outputScoreMtx(ScoreMtx);	
	//cout<< "finish initialization"<<endl;
	reset_cor(Coors);
	vector<int> nzero;	
	RecCal(0,Coors,ScoreMtx,nzero);	
	//cout<<"finish calculation"<<endl;
	//outputScoreMtx(ScoreMtx);
	TraceAln(ScoreMtx);	
	//cout<<"Length is "<<seq1.size()
	this->FinalScore= ScoreMtx[TotalSize-1].val/(float)this->otseqs[0].size();
	delete[] ScoreMtx;
	delete[] Coors;
	return this->FinalScore;
}

//trace back and form the actual alignment
void Alignment::TraceAln(MtxE* mtx){
	for(unsigned int i=0; i< this->inseqs.size();i++){
		vector<string> temp;
		this->otseqs.push_back(temp);
	}
	int* Cors= new int[(int)this->inseqs.size()];
	for(unsigned int i=0; i < this->inseqs.size();i++)
		Cors[(int)i]=(int)this->inseqs[i].size();

	
	int index_Cors = cor_to_ind(Cors);
	int shift;
	while(mtx[index_Cors].prev > 0){
		for(unsigned int i=0; i < this->inseqs.size();i++){
			shift = (int)this->inseqs.size()-i-1;
			if((mtx[index_Cors].prev & (1<<shift)) == 0 )
				this->otseqs[i].insert(this->otseqs[i].begin(),"-");
			else{
				this->otseqs[i].insert(this->otseqs[i].begin(),this->inseqs[i][Cors[i]-1]);
				Cors[i]--;
				
			}
		}
		index_Cors = cor_to_ind(Cors);	
	}	
	delete[] Cors;		
	return;
}

void Alignment::RecCal(int lev,int* Cors,MtxE* mtx,vector<int> nzero){
	if(lev < (int)this->inseqs.size()){
		Cors[lev]=0;
		RecCal(lev+1,Cors,mtx,nzero);
		nzero.push_back(lev);
		for(int i=1; i <= (int)this->inseqs[lev].size();i++){
			Cors[lev]=i;
			RecCal(lev+1,Cors,mtx,nzero);
		}

	}else{
		if((int)nzero.size() < 2) return;//base cases are already taken care
		int* NCors = new int[(int)this->inseqs.size()];
		//use a int to encode all possible neighbors
		int allnbr= 1<<(int)(this->inseqs.size());
		int shift;//use to store the number of shift needed for bitwise ops
		int mask;//mask used to check the validity of neighbor positions
		bool first=1;//1 if the score of current position has not set yet.
		//to record the alignment type at Cors position for each seq,
		//1 means it is a gap, 0 otherwise
		int index_Cors=cor_to_ind(Cors);
		int index_NCors;
		bool* alntyp = new bool[(int)this->inseqs.size()];
		for(int cur=1; cur < allnbr;cur++){
			for(int j=0; j < (int)this->inseqs.size();j++)
				NCors[j] = Cors[j];
			mask = 0;
			for(unsigned int j=0;j < nzero.size();j++){
				shift = (int)this->inseqs.size()-nzero[j]-1;
				mask+=(1<<shift);
			}

			/**FOR DEBUG see the bit string of the number*/
			//std::bitset<3> temp2(mask);
			//cout<<"M:"<<temp2<<endl;
			//std::bitset<3> temp1(cur);
			//cout<<"N:"<<temp1<<endl;

			//some neighbor positions may not exist
			//for some element, skip cur neighbor if not exit.
			if( (cur|mask) != mask) continue;

			for(int j=0; j < (int)this->inseqs.size();j++){
				shift = (int)this->inseqs.size()-j-1;
				if((cur & (1<<shift)) > 0){
					 NCors[j]--;
					alntyp[j]=0;
				}else{
					alntyp[j]=1;
				}
			}
			/****For Debug**********
			cout<<"---"<<endl;	
			for(int j=0; j < (int)this->inseqs.size();j++)
				cout<<Cors[j]<<' ';
			cout<<endl;
			for(int j=0; j < (int)this->inseqs.size();j++)
				cout<<NCors[j]<<' ';
			cout<<endl;
			for(int j=0; j < (int)this->inseqs.size();j++)
				cout<<alntyp[j]<<' ';
			cout<<endl;
			***********************/

			index_NCors = cor_to_ind(NCors);
			
			if(mtx[index_NCors].prev == -1)
				cout<<"Error!!!this NCors is not computed yet!!"<<endl;
	
			float cur_score = mtx[index_NCors].val+pos_score(mtx,Cors,alntyp);
			//cout<<cur_score<<endl;	
			if(first == 1){
				mtx[index_Cors].val = cur_score;
				mtx[index_Cors].prev = cur;
				first = 0;
			}else if(mtx[index_Cors].val > cur_score){
				mtx[index_Cors].val = cur_score;
				mtx[index_Cors].prev = cur;
			}	

		}
		delete[] NCors;
		delete[] alntyp;
	}
	return;
}
//1 means it is a gap, 0 otherwise
//assume the length of Cors, Type must be of length this->inseqs.size()
float Alignment::pos_score(MtxE* mtx,int* Cors,bool* Type){
	float score_sum=0;
	for(unsigned int i=0; i < this->inseqs.size();i++){
		if(Type[i] == 1){
			score_sum += this->Gap*((int)this->inseqs.size()-i-1);
			 continue;
		}
		for(unsigned int j=i+1;j < this->inseqs.size();j++){
			if(Type[j] == 1)
				score_sum += this->Gap;
			else{		
				float Percentage=Similarity(this->inseqs[i][Cors[i]-1],this->inseqs[j][Cors[j]-1]);
				if(FloatEqual(Percentage,1.0) == true ){
					score_sum+=this->Match;
				}else if(FloatEqual(Percentage,0.9) == true){
					score_sum+=this->Identical;
				}else if(Percentage >= 0.66 ){
					score_sum+=this->Similar;
				}else{
					score_sum+=this->Mismatch;
				}
			}
		}
	}
	return (2.0*score_sum)/(((int)this->inseqs.size())*((int)this->inseqs.size()-1));
}
	
void Alignment::outputScoreMtx(MtxE* mtx){
	int* Coors= new int[(int)this->inseqs.size()];
	reset_cor(Coors);
	cout<<"\t0";
	for(unsigned int i=0; i < this->inseqs.back().size();i++)
		cout<<"\t"<<this->inseqs.back()[i];
	cout<<endl;
	outputScoreMtxHelper(0,Coors,mtx);
	delete[] Coors;			
	return;
}

void Alignment::outputScoreMtxHelper(int lev,int* Cors,MtxE* mtx){
	if(lev < (int)this->inseqs.size()-1){
		for(int i=0; i <= (int)this->inseqs[lev].size();i++){
			Cors[lev]=i;
			outputScoreMtxHelper(lev+1,Cors,mtx);
		}
	}else{
		if(Cors[0] == 0)
			cout<<"0";
		else
			cout<<this->inseqs[0][Cors[0]-1];
		for(int i=1; i < (int)this->inseqs.size()-1;i++)
			if(Cors[i] == 0)
				cout<<"0";
			else
				cout<<","<<this->inseqs[i][Cors[i]-1];

		for(int i=0; i <= (int)this->inseqs[lev].size();i++){
			Cors[lev]=i;
			cout<<"\t"<<mtx[cor_to_ind(Cors)].val<<"/"<<mtx[cor_to_ind(Cors)].prev;
		}
		cout<<endl;
	}
	return;
}
	
void Alignment::PrintAlignment(){
	cout<<this->FinalScore<<endl;
	for(unsigned int i=0; i < this->otseqs.size();i++){
		cout<<this->names[i];
		for(unsigned int j=0; j < this->otseqs[i].size();j++)
			cout<<"\t"<<this->otseqs[i][j];
		cout<<endl;
	}
	return;
}
void Alignment::PrintAlignment(vector<string> & store, int index){
//	if(index == 1){
//		for(unsigned int j=0; j <seq1a.size();j++){
//			store.push_back(seq1a[j]);
//		}
//	}else if(index == 2){
//		for(unsigned int j=0; j <seq2a.size();j++){
//			store.push_back(seq2a[j]);
//		}
//	}
	return;
}
		
	

