#include "align.h"

Alignment::Alignment(vector<string>& a,vector<string>& b){
	name1=a[0];
	name2=b[0];
	for(unsigned int i=1;i < a.size();i++){
		seq1.push_back(a[i]);
	}
	for(unsigned int i=1;i < b.size();i++){
		seq2.push_back(b[i]);
	}
	gap=1.0;
	mismatch=1.0;
	identical=-1.0;
	similar= -0.5;
	match= -2.0;
}

Alignment::Alignment(vector<string>& a,string namea,vector<string>& b,string nameb){
	name1=namea;
	name2=nameb;
	for(unsigned int i=0;i < a.size();i++){
		seq1.push_back(a[i]);
	}
	for(unsigned int i=0;i < b.size();i++){
		seq2.push_back(b[i]);
	}
	//cout<<"Created a alignment of size "<<seq1.size()<<" AND "<<seq2.size()<<endl;
	gap=1.0;
	mismatch=1.0;
	identical=-1.0;
	similar= -0.5;
	match= -2.0;
}
	

Alignment::Alignment(vector<string>& a,vector<string>& b, float c,float d,float e,float f,float g){
	name1=a[0];
	name2=b[0];
	for(unsigned int i=1;i < a.size();i++){
		seq1.push_back(a[i]);
	}
	for(unsigned int i=1;i < b.size();i++){
		seq2.push_back(b[i]);
	}
	gap=c;
	mismatch=d;
	identical=e;
	similar=f;
	match=g;
}

Alignment::~Alignment(){
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


float Alignment::align(void){
	float  temp,min;
	int ydim=seq1.size()+1;
	int xdim=seq2.size()+1;
	float* score = new float[xdim*ydim];
	for(int i=0; i < ydim;i++){
		score[i*xdim]=(float)i*gap;
	}
	for(int i=0; i <xdim;i++){
		score[i]=(float)i*gap;
	}
	
	//cout<< "finish initialization"<<endl;
	
	for(int j=1; j <xdim;j++){
		for(int i=1; i <ydim;i++){
			min=score[(i-1)*xdim+j]+gap;
		
			//seq1a[m]=seq1[m];
			//seq2a[n]="-";
			float Percentage=Similarity(seq1[i-1],seq2[j-1]);
			
			if(FloatEqual(Percentage,1.0) == true ){
				temp=score[(i-1)*xdim+j-1]+match;
			}else if(FloatEqual(Percentage,0.9) == true){
				temp=score[(i-1)*xdim+j-1]+identical;
				//cout<<"Identical "<<temp<<endl;
			}else if(Percentage >= 0.66 ){
				temp=score[(i-1)*xdim+j-1]+similar;
				//cout<<"Similar "<<temp<<endl;
			}else{
				temp=score[(i-1)*xdim+j-1]+mismatch;
				//cout<<"Mismatch "<<temp<<endl;
			}
			if(temp < min){
				min=temp;
			}
			temp=score[i*xdim+j-1]+gap;
			if(temp < min){
				min=temp;
			}
			//cout<<seq1[i-1]<<","<<seq2[j-1]<<":"<<Percentage<<",Score:"<<min<<endl;
			score[i*xdim+j]=min;
		}
	}
	//cout<<"finish calculation"<<endl;
	
	/**************  Output score matrix*************
	for(int i=ydim-1; i >=0;i--){
		if( i >= 1){
			cout<<seq1[i-1];
		}else{
			cout<<"-";
		}
		for(int j=0; j <xdim;j++){
			cout<<"\t"<<score[i*xdim+j];
		}
		cout<<endl;
	}
	cout<<"\t"<<"-";
	for(unsigned int j=0; j < seq2.size();j++){
		cout<<"\t"<<seq2[j];
	}
	cout<<endl;
	cout<<"Final Score :"<<score[xdim*ydim-1]<<endl;
	*************************************************/
	int xInd=xdim-1;
	int yInd=ydim-1;
	while(xInd >=1 || yInd >=1 ){
		//cout <<"Checking "<<xInd<<","<<yInd<<endl;
		if(xInd-1 >=0 && score[yInd*xdim+xInd] == score[yInd*xdim+xInd-1]+gap){
			seq1a.insert(seq1a.begin(),"-");
			seq2a.insert(seq2a.begin(),seq2[xInd-1]);
			xInd=xInd-1;
			continue;
		}
		if(yInd-1 >=0 && score[yInd*xdim+xInd] == score[(yInd-1)*xdim+xInd]+gap){
			seq1a.insert(seq1a.begin(),seq1[yInd-1]);
			seq2a.insert(seq2a.begin(),"-");
			yInd=yInd-1;
			continue;
		}
		if(score[yInd*xdim+xInd] == score[(yInd-1)*xdim+xInd-1]+mismatch || score[yInd*xdim+xInd] == score[(yInd-1)*xdim+xInd-1]+identical || score[yInd*xdim+xInd] == score[(yInd-1)*xdim+xInd-1]+similar || score[yInd*xdim+xInd] == score[(yInd-1)*xdim+xInd-1]+match){
			seq1a.insert(seq1a.begin(),seq1[yInd-1]);
			seq2a.insert(seq2a.begin(),seq2[xInd-1]);
			xInd=xInd-1;
			yInd=yInd-1;
		}
	}
	//cout<<"Length is "<<seq1.size()
	float n_score = score[xdim*ydim-1]/(float)seq1a.size();
	delete[] score;
	score = NULL;
	return n_score;
}

void	 Alignment::outputalign(){
	cout<<name1;
	for(unsigned int j=0; j <seq1a.size();j++){
		cout<<"\t"<<seq1a[j];
	}
	cout<<endl;
	cout<<name2;
	for(unsigned int j=0; j <seq2a.size();j++){
		cout<<"\t"<<seq2a[j];
	}
	cout<<endl;
	return;
}
void Alignment::outputalign(vector<string> & store, int index){
	if(index == 1){
		for(unsigned int j=0; j <seq1a.size();j++){
			store.push_back(seq1a[j]);
		}
	}else if(index == 2){
		for(unsigned int j=0; j <seq2a.size();j++){
			store.push_back(seq2a[j]);
		}
	}
	return;
}
		
	
//~ float Alignment::subalign(int m,int n){
	//~ float min,temp;
	//~ if(m==0 || n==0){
		//~ return score[m*seq2.size()+n];
	//~ }
	//~ min=subalign(m-1,n)+gap;
		
			//~ //seq1a[m]=seq1[m];
			//~ //seq2a[n]="-";
	//~ if(seq1[m].compare(seq2[n]) == 0){
		//~ temp=subalign(m-1,n-1);
	//~ }else{
		//~ temp=subalign(m-1,n-1)+mismatch;
	//~ }
			
	//~ if(temp < min){
		//~ min=temp;
	//~ }
	//~ temp=subalign(m,n-1)+gap;
	//~ if(temp < min){
		//~ min=temp;
	//~ }
			
			
	//~ //cout<<m<<";"<<n<<endl;
	//~ //cout<<gap<<";"<<mismatch<<endl;
	//~ return min;
//~ }
//		void Alignment::outputalign();
