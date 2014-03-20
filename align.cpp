#include "align.h"

Alignment::Alignment(){	
	Gap=1.0;
	Mismatch=1.0;
	Identical=-1.0;
	Similar= -0.5;
	Match= -2.0;
}

Alignment::Alignment(float c,float d,float e,float f,float g){
	Gap=c;
	Mismatch=d;
	Identical=e;
	Similar=f;
	Match=g;
}

Alignment::~Alignment(){
}

bool Alignment::AddSeq(vector<string> const & input){
	if(input.size() < 2){
		cerr<<"Seq to add too short, it needs to contain at least a name and a single string "<<endl;
		return 0;
	}
	vector<string> newseq = input;
	this->names.push_back(newseq.front());
	newseq.erase(newseq.begin());
	if(this->inseqs.size() == 0 || newseq.size() >= this->inseqs[0].size()  )
		this->inseqs.push_back(newseq);
	else
		this->inseqs.insert(this->inseqs.begin(),newseq);
	
	return 1;
}

bool Alignment::AddSeq(vector<string> const & input,string const & iname){
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
	if(this->inseqs.size() == 0 || newseq.size() >= this->inseqs[0].size()  )
		this->inseqs.push_back(newseq);
	else
		this->inseqs.insert(this->inseqs.begin(),newseq);
	
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
	for(int i=0;i < this->MaxSeq;i++){
		val[i]=0;
	}	
}
		
float Alignment::Align(void){
	float  temp,min;
	if((int)this->inseqs.size() < 2){
		cerr<<"There is only 1 seq instored,can't do alignment"<<endl;
		return 10;
	}
	int TotalSize=1;
	for(unsigned int i=0;i < this->inseqs.size();i++)
		TotalSize = TotalSize*((int)this->inseqs[i].size()+1);

	MtxE * ScoreMtx = new MtxE[TotalSize];
	int Coors[this->MaxSeq];
	reset_cor(Coors);
	
	for(unsigned int i=0; i < this->inseqs.size();i++){
		for(int j=0; j <= (int)this->inseqs[i].size();j++){
			if(i > 0) Coors[i-1]=0;
			Coors[i]=j;
			ScoreMtx[cor_to_ind(Coors)].val=(float)(j*this->Gap);
			ScoreMtx[cor_to_ind(Coors)].prev = -2;	
		}
	}	
	outputScoreMtx(ScoreMtx);	
	//cout<< "finish initialization"<<endl;
	reset_cor(Coors);
	vector<int> nzero;	
	RecCal(0,Coors,ScoreMtx,nzero,2);	
	//cout<<"finish calculation"<<endl;
		
	//cout<<"Length is "<<seq1.size()
	float n_score = 0.1;//ScoreMtx[TotalSize-1].val/(float)this->otseqs[0].size();
	delete[] ScoreMtx;
	ScoreMtx = NULL;
	return n_score;
}


void Alignment::RecCal(int lev,int* Cors,MtxE* mtx,vector<int> nzero,int limit){
	if(lev < (int)this->inseqs.size() && (int)nzero.size() <= limit){
		Cors[lev]=0;
		RecCal(lev+1,Cors,mtx,nzero,limit);
		if((int)nzero.size() < limit){
			nzero.push_back(lev);
			for(int i=1; i < (int)this->inseqs[lev].size();i++){
				Cors[lev]=i;
				RecCal(lev+1,Cors,mtx,nzero,limit);
			}
		}
	}else{
		if((int)nzero.size() < 2) return;//base cases are already taken care
		int NCors[this->MaxSeq];
		float min_score;
		reset_cor(NCors);
		for(int i=0; i < )
		//CalScore(NCors,Cors,)
		
	}
	return;
}
	
void Alignment::outputScoreMtx(MtxE* mtx){	
	int Coors[this->MaxSeq]={0};
	cout<<"\t0";
	for(unsigned int i=0; i < this->inseqs.back().size();i++)
		cout<<"\t"<<this->inseqs.back()[i];
	cout<<endl;
	outputScoreMtxHelper(0,Coors,mtx);	
	return;
}

void Alignment::outputScoreMtxHelper(int lev,int* Cors,MtxE* mtx){
	if(lev < (int)this->inseqs.size()-1){
		for(int i=0; i < (int)this->inseqs[lev].size();i++){
			Cors[lev]=i;
			outputScoreMtxHelper(lev+1,Cors,mtx);
		}
	}else{
		cout<<Cors[0];
		for(int i=1; i < (int)this->inseqs.size()-1;i++)
			cout<<","<<Cors[i];

		for(int i=0; i < (int)this->inseqs[lev].size();i++){
			Cors[lev]=i;
			cout<<"\t"<<mtx[cor_to_ind(Cors)].val;
		}
		cout<<endl;
	}
	return;
}
	
void Alignment::PrintAlignment(){
//	cout<<name1;
//	for(unsigned int j=0; j <seq1a.size();j++){
//		cout<<"\t"<<seq1a[j];
//	}
//	cout<<endl;
//	cout<<name2;
//	for(unsigned int j=0; j <seq2a.size();j++){
//		cout<<"\t"<<seq2a[j];
//	}
//	cout<<endl;
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
		
	

