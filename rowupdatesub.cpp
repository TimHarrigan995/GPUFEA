//============================================================================
// Name        : rowupdate.cpp
// Author      : Tim
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <vector>
#include <array>
#include <bits/stdc++.h>
#include <numeric> // For std::accumulate (optional, for calculating total_size)


using namespace std;

// Function for inserting elements
// in array of vectors
void InsertionIntoArray(vector <int> &vdofs, int doftoinsert,int row)
{
int c;
		if(doftoinsert<row)return;
		if(vdofs.size()==0){
            vdofs.push_back(doftoinsert);
            return;
		}
		else {

			for (c=0;c<vdofs.size();c++){
				if(doftoinsert == vdofs[c])return;
				if(doftoinsert < vdofs[c]){
					vdofs.insert(vdofs.begin()+c,doftoinsert);
					return;
				} // if(dof...
			} //for..
			if(doftoinsert > vdofs[vdofs.size()-1])vdofs.push_back(doftoinsert);
			//getting to this statement means doftoinsert is

		}
            return;
} //

void elementDOFconnectivity(int elgdofs[], int elsiz,vector <int>* globaldofs,int aSize){

	//set up calls to insert array dofs in the row arrays
	for(int growind=0;growind<elsiz;growind++){
		int globalrow=elgdofs[growind]-1;

		if(globalrow>=0){
		for (int c=0;c<elsiz;c++){
			if(elgdofs[c]>0){
			InsertionIntoArray(globaldofs[globalrow],elgdofs[c]-1,globalrow);}
		} //columns (=dof-1)
		}
	} //rows (=dof-1)

}

/*
int main() {


    // Declare a pointer to an array of row vectors std::vector<int>
    std::vector<int>* gldofs;
    std::vector<int> StrungDofs;

	gldofs = new std::vector<int>[100];

	int arrSize=100;
	int diagindex[100];

	int elgsiz=8;
	int elgdofs[8]={3,4,5,6,13,14,11,12};

	elementDOFconnectivity(elgdofs, elgsiz, gldofs,arrSize);

	int elgdofs2[8]={1,2,3,4,11,12,9,10};

	elementDOFconnectivity(elgdofs2, elgsiz, gldofs,arrSize);

	int elgdofs3[8]={0,0,0,0,3,4,1,2};
	elementDOFconnectivity(elgdofs3, elgsiz, gldofs,arrSize);

	int elgdofs4[8]={0,0,0,0,5,6,3,4};
	elementDOFconnectivity(elgdofs4, elgsiz, gldofs,arrSize);

	int elgdofs5[8]={5,6,7,8,15,16,13,14};
	elementDOFconnectivity(elgdofs5, elgsiz, gldofs,arrSize);


//concatenating vectors to a single full vector

	size_t total_size = 0;
    for (int c=0;c<17;c++) {
    	diagindex[c]=total_size;
        total_size += gldofs[c].size();

    std::cout << "Vector elements: ";
    for (int num : gldofs[c]) {
        std::cout << num << " ";
    }
    std::cout << std::endl;
    }



    StrungDofs.reserve(total_size); // Pre-allocate memory for efficiency

    for (int c=0;c<17;c++) {
        StrungDofs.insert(StrungDofs.end(), gldofs[c].begin(), gldofs[c].end());}

    std::cout << "Vector elements: ";
    for (int num : StrungDofs) {
        std::cout << num << " ";
    }
    std::cout << std::endl;


    return 0;

}
*/
