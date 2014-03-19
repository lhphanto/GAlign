#include "malign.h"

MAlignment::MAlignment(float a,float b,float c){
	gap=a;
	mismatch=b;
	similar=c;
}
void MAlignment::AddSeq(vector<string>& a){
	names.push_back(a[0]);	
	vector<string> temp;
	for(unsigned int i=1;i < a.size();i++){
		temp.push_back(a[i]);
	}
	seqs.push_back(temp);
	return;
}
void MAlignment::AddSeq(vector<vector<string> > &  a){
	for(unsigned int i=0; i < a.size();i++){
		AddSeq(a[i]);
	}
	return;
}


float MAlignment::align(){
	profile.push_back(seqs[0]);
	float FinalScore;
	int total_length=seqs[0].size();
	
	for(unsigned int i=1;i < seqs.size();i++){
		total_length+=seqs[i].size();
		/*******Calculate the TABLE********/		
		int ydim=profile[0].size()+1;
		int xdim=seqs[i].size()+1;
		float score[xdim*ydim];
		int trace[xdim*ydim]; // to record traces for finding the aligned sequences 1 for gap in x,2 for match,3 for gap in y
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
				score[k*xdim+j]=score[(k-1)*xdim+j]+gap;//initialize min with adding gap in y direction
				trace[k*xdim+j]=3;
				for(unsigned int m=0;m < profile.size();m++){
					float Percentage=Similarity(profile[m][k-1],seqs[i][j-1]);
					if(FloatEqual(Percentage,1.0) == true){
						temp=temp+score[(k-1)*xdim+j-1]+identical;
						//cout<<"Identical "<<profile[m][k-1]<<","<<seqs[i][j-1]<<" Score:"<<temp<<endl;
					}else if(FloatEqual(Percentage,0.9) == true){						
						temp=temp+score[(k-1)*xdim+j-1]+verysimilar;
						//cout<<"Almost "<<profile[m][k-1]<<","<<seqs[i][j-1]<<" Score:"<<temp<<endl;
					}else if(Percentage >= 0.66 ){
						//cout<<"Temp now is "<<temp<<endl;
						temp=temp+score[(k-1)*xdim+j-1]+similar;
						//cout<<"Similar "<<profile[m][k-1]<<","<<seqs[i][j-1]<<" Score:"<<temp<<endl;
					}else if(Percentage < 0.0 ){
						temp=temp+score[(k-1)*xdim+j-1]+gap;
						//cout<<"MatchToGap "<<profile[m][k-1]<<","<<seqs[i][j-1]<<" Score:"<<temp<<endl;
					}else{						
						temp=temp+score[(k-1)*xdim+j-1]+mismatch;
						//cout<<"Mismatch "<<profile[m][k-1]<<","<<seqs[i][j-1]<<" Score:"<<temp<<endl;
					}
				}
				temp=temp/float(profile.size());
				//cout<<"Finally the score is "<<temp<<endl;
				if(temp < score[k*xdim+j]){
					//cout<<temp << " is smaller than "<< score[k*xdim+j] <<endl;
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
		while(xInd > 0 || yInd > 0 ){
			//cout <<"Checking "<<xInd<<","<<yInd<<endl;
			if(xInd-1 >=0 && trace[yInd*xdim+xInd] == 1){
				for(unsigned int m=0;m < profile.size();m++){
					if(alignment.size() > m){
						alignment[m].insert(alignment[m].begin(),"-");
					}else{
						vector<string> temp;
						temp.push_back("-");
						alignment.push_back(temp);
					}
				}
					
				if(alignment.size() > profile.size() ){
					alignment[int(profile.size())].insert(alignment[int(profile.size())].begin(),seqs[i][xInd-1]);
				}else{
					vector<string> temp;
					temp.push_back(seqs[i][xInd-1]);
					alignment.push_back(temp);
				}
				xInd=xInd-1;
				continue;
			}
			
			if(yInd-1 >=0 && trace[yInd*xdim+xInd] == 3){
				for(unsigned int m=0;m < profile.size();m++){
					if(alignment.size() > m){
						alignment[m].insert(alignment[m].begin(),profile[m][yInd-1]);
					}else{
						vector<string> temp;
						temp.push_back(profile[m][yInd-1]);
						alignment.push_back(temp);
					}
				}
					
				if(alignment.size() > profile.size()){
					alignment[int(profile.size())].insert(alignment[int(profile.size())].begin(),"-");
				}else{
					vector<string> temp;
					temp.push_back("-");
					alignment.push_back(temp);
				}			
				yInd=yInd-1;
				continue;
			}
			if(trace[yInd*xdim+xInd] == 2){
				for(unsigned int m=0;m < profile.size();m++){
					if(alignment.size() > m){
						alignment[m].insert(alignment[m].begin(),profile[m][yInd-1]);
					}else{
						vector<string> temp;
						temp.push_back(profile[m][yInd-1]);
						alignment.push_back(temp);
					}
				}
					
				if(alignment.size() > profile.size()){
					alignment[int(profile.size())].insert(alignment[int(profile.size())].begin(),seqs[i][xInd-1]);
				}else{
					vector<string> temp;
					temp.push_back(seqs[i][xInd-1]);
					alignment.push_back(temp);
				}				
				xInd=xInd-1;
				yInd=yInd-1;
			}
		}
		
		 //~ while (i > 0)
  //~ {
    //~ AlignmentA <- Ai + AlignmentA
    //~ AlignmentB <- "-" + AlignmentB
    //~ i <- i - 1
  //~ }
  //~ while (j > 0)
  //~ {
    //~ AlignmentA <- "-" + AlignmentA
    //~ AlignmentB <- Bj + AlignmentB
    //~ j <- j - 1
  //~ }

		UpdateProfile();
		FinalScore=score[xdim*ydim-1];
	}
	
	return (FinalScore)/(int)(profile[0].size())+2;
	
}

void MAlignment::UpdateProfile(){
	for(unsigned int i=0; i< profile.size();i++){
		profile[i].clear();
		for(unsigned int j=0;j < alignment[i].size();j++){
			string temp = alignment[i][j];
			profile[i].push_back(temp);
		}
	}
	profile.push_back(alignment[int(profile.size())]);
	for(unsigned int j=0;j < alignment.size();j++){
		alignment[j].clear();
	}
	return;
}
	
void  MAlignment::outputalign(){
	for(unsigned int i=0; i< profile.size();i++){
		cout<<names[i];
		for(unsigned int j=0;j < profile[i].size();j++){
			cout<<"\t"<<profile[i][j];
		}
		cout<<endl;
	}
	return;
}

//Helper function to compare two float numbers to 0.01 precision
bool MAlignment::FloatEqual(float A, float B){
    if (A == B)
        return true;

    float Error = fabs((A - B));
    if (Error <= 0.01)
        return true;
    return false;

}
		
float  MAlignment::Similarity(string a,string b){
	string tempa,tempb;
	float sim=0.0;
	if(a.size() != b.size() || a.compare("-") == 0 || b.compare("-") == 0){
		return -1.0;
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

	if( sim == ((float)a.size()-0.5) ){
		return 0.9;
	}

	if( sim == (float)a.size() ){
		return 1;
	}
	
		
	return ((float)((int)sim))/(float)a.size();	
}
