/**
 * Shape Cut for Conjugate Product Graph
 *
 */

#include "mex.h"

#include <exception>
#include <fstream>
#include <list>
#include <iostream>
#include <string>
#include <cmath>

#include "buildMeshCGapprox.hpp"
#include "helperCGapprox.hpp"
#include "manifoldDijkstraCGapprox.hpp"
#include "MinHeapCGapprox.hpp"
#include "shapeCutCGapprox.hpp"
#include "dijkstraCGapprox_main.hpp"

using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
// check and read mex input
	if (nrhs != 6 || nlhs != 3)
		mexErrMsgTxt("Usage: [energy, path, correspondence] = dijkstra(vertices, triangles, contour, landmarks_2d, landmarks_3d, quaternions). Transform indices in triangles to initial 0.");

	// read vertices
	const double* const vertices = mxGetPr(prhs[0]);
	const long nVertices = long( mxGetN(prhs[0]) );
	const int dimVertices = int( mxGetM(prhs[0]) );

	// read triangles
	const double* const triangles = mxGetPr(prhs[1]);
	const long nTriangles = long( mxGetN(prhs[1]) );
	const int dimTriangles = int( mxGetM(prhs[1]) );

	// read contour
	const double* const contour = mxGetPr(prhs[2]);
	const long nContour = long( mxGetN(prhs[2]) ) - 1; // because we duplicated two layers instead of one (same in CG code)
	const int dimContour = int( mxGetM(prhs[2]) );
    
    int32_T* landmarks_2d = (int32_T *)mxGetPr(prhs[3]);
	const int32_T* landmarks_3d = (int32_T *)mxGetPr(prhs[4]);
    int nLandmarks = int( mxGetN(prhs[3]) );

    if(int( mxGetN(prhs[3]) )!=int( mxGetN(prhs[4]) )) nLandmarks = -1;


    double *matrixData = mxGetPr(prhs[5]);
    size_t numQuaternions = mxGetN(prhs[5]);

    std::vector<Eigen::Quaterniond> precRotations;
    
    if(numQuaternions!=0)
    {
        precRotations.reserve(numQuaternions);
        
        for (size_t i = 0; i < numQuaternions; ++i) {
            double w = matrixData[i*4 + 0];
            double x = matrixData[i*4 + 1];
            double y = matrixData[i*4 + 2];
            double z = matrixData[i*4 + 3];
            precRotations.emplace_back(w, x, y, z);
        }
    }
    
    // cost mode 
    const double costMode = 13;

	if (dimTriangles != 3)
		mexErrMsgTxt("Triangles must be given in 3xm.");

	//if (dimContour != dimVertices)
	//	mexErrMsgTxt("Descriptors must have the same dimension. (Remove this constraint if you changed the energy function to handle different dimensional features.)");

try
{

    std::vector<std::tuple<long,long>> matching;
	const double energy = dijkstra_main(vertices, nVertices, dimVertices, triangles, nTriangles, contour, nContour, dimContour, matching, landmarks_2d, landmarks_3d, nLandmarks, precRotations);
    
	// transform results for mex

	// return loop
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	plhs[1] = mxCreateDoubleMatrix(matching.size(), 1, mxREAL);
	plhs[2] = mxCreateDoubleMatrix(matching.size(), 1, mxREAL);
	double* e = mxGetPr(plhs[0]);
	double* resultPath = mxGetPr(plhs[1]);
	double* resultCorr = mxGetPr(plhs[2]);

    if(nLandmarks == -1){
        return;
    }
	// energy
	e[0] = energy;

	// path
    int i = 0;
	for (const auto& match : matching) {
		resultPath[i] = std::get<0>(match); // 3d vertex id
		resultCorr[i] = std::get<1>(match); // 2d vertex id
        i++;
	}
    resultCorr[i-1]--;

} catch (std::exception& e)
	{
		const std::string msg = "Exception caught: " + std::string(e.what());
		mexErrMsgTxt(msg.c_str());
	}

} // close mexFunction
