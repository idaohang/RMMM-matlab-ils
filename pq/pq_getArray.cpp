//==============================================================================
// Name        : pq_top.cpp
// Author      : Stephen Breen
// Version     : 1.0
// Description : Returns all the elements currently in the heap as an array
//
// May 22, 2009: Created
//==============================================================================
#include "MyHeap.h"

//------------------------------- MATLAB -------------------------------------//
 #define toSysout(...) printf(__VA_ARGS__)
 #define exit_with_error(...)           \
 do {                                   \
		 fprintf(stdout, "Error: ");    \
		 fprintf(stdout, __VA_ARGS__ ); \
		 fprintf(stdout, "\n" );        \
		 exit(1);                       \
 } while(0)
#ifdef MATLAB_MEX_FILE
#include "mex.h"

void retrieve_heap( const mxArray* matptr, MinHeap<double>* & heap){
    // retrieve pointer from the MX form
    double* pointer0 = mxGetPr(matptr);
    // check that I actually received something
    if( pointer0 == NULL )
        mexErrMsgTxt("vararg{1} must be a valid priority queue pointer\n");
    // convert it to "long" datatype (good for addresses)
    long pointer1 = (long) pointer0[0];
    // convert it to "KDTree"
    heap = (MinHeap<double>*) pointer1;
    // check that I actually received something
    if( heap == NULL )
        mexErrMsgTxt("vararg{1} must be a valid priority queue pointer\n");
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if( nrhs!=1 )
		mexErrMsgTxt("This function requires 1 arguments\n");
	if( !mxIsNumeric(prhs[0]) )
		mexErrMsgTxt("parameter 1 missing!\n");

	// retrieve the heap
	MinHeap<double>*  heap;
	retrieve_heap( prhs[0], heap);
	vector<pair<double, int> > * heapArr = heap->getArray();
    //Allocate memory and assign output pointer
    plhs[0] = mxCreateDoubleMatrix(1, heap->size(), mxREAL); //mxReal is our data-type
    //Get a pointer to the data space in our newly allocated memory
    double *outArray = mxGetPr(plhs[0]);
    int i;
    for(i=0;i<heapArr->size();i++)
    {
        outArray[i] = heapArr->at(i).second+1;
    }
}
#endif
