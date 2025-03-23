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

#include "buildMeshCG.hpp"
#include "helperCG.hpp"
#include "arap.hpp"


using namespace std;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
// check and read mex input
	if (nrhs != 3 || nlhs != 1)
		mexErrMsgTxt("Usage: [rotations] = dijkstra(vertices, triangles, contour). Transform indices in triangles to initial 0.");

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
    

    
    // cost mode 
    const double costMode = 13;

	if (dimTriangles != 3)
		mexErrMsgTxt("Triangles must be given in 3xm.");

	//if (dimContour != dimVertices)
	//	mexErrMsgTxt("Descriptors must have the same dimension. (Remove this constraint if you changed the energy function to handle different dimensional features.)");

try
{
    vector<list<long>* > mesh(nVertices);
    vector<list<std::tuple<long,long>>* > mesh2(nVertices);
	for (long i=0; i<nVertices; i++) {
		mesh[i] = new list<long>();
        mesh2[i] = new list<std::tuple<long, long>>();
	}
    vector<std::tuple<long,long>> edges3d;
	buildMesh(triangles, nTriangles, mesh, mesh2, edges3d);

    
    const long nEdges3d = edges3d.size();
    const long nProductNodesPerLayer = 6 * nEdges3d;
    const long nProductNodes = nProductNodesPerLayer * 3;

    

    LGProdNode nodesPerLayer[nProductNodesPerLayer];
    for (int i = 0; i < nProductNodesPerLayer; i++) {
        nodesPerLayer[i] = getProductNodeInLayer(i, 0, edges3d, nVertices);
    }
    std::vector<Eigen::Quaterniond> precRotations;

    if (PRECOMPUTE_ROTATIONS){
        precRotations.reserve(nProductNodesPerLayer*(nContour+1));
        for (int i = 0; i < nContour+1; i++) {
            for (int j = 0; j < nProductNodesPerLayer; j++) {
                const long idx = nContour * nProductNodesPerLayer + j;
                LGProdNode cn = nodesPerLayer[j];
                // add level
                cn.idxSrc2d += i; cn.idxTrgt2d += i;
                long c2dTrgt = cn.idxTrgt2d;
                long c2dSrc = cn.idxSrc2d;
                c2dSrc =  cn.idxSrc2d == cn.idxTrgt2d  ? cn.idxSrc2d-1 : cn.idxSrc2d;
                if (c2dSrc < 0) c2dSrc = nVertices-2;
                if (c2dSrc == c2dTrgt) {
                    std::vector<int> throwVec;
                    throwVec.at(-1);
                }
                const Eigen::Quaterniond rot = precomputeRotation(vertices, nVertices,dimVertices,
                        contour,nContour,dimContour,c2dSrc, c2dTrgt, cn.idxSrc3d, cn.idxTrgt3d);
                precRotations.push_back(rot);
            }
        }
    }

    
	// transform results for mex

    mxArray *quaternionMatrix = mxCreateDoubleMatrix(4, nProductNodesPerLayer*(nContour+1), mxREAL);
    double *matrixData = mxGetPr(quaternionMatrix);
    
    // Populate the matrix with the quaternion data
    for (size_t i = 0; i < nProductNodesPerLayer*(nContour+1); ++i) {
        matrixData[i*4 + 0] = precRotations[i].w();
        matrixData[i*4 + 1] = precRotations[i].x();
        matrixData[i*4 + 2] = precRotations[i].y();
        matrixData[i*4 + 3] = precRotations[i].z();
    }
    
    // Return the matrix to MATLAB
    plhs[0] = quaternionMatrix;

} catch (std::exception& e)
	{
		const std::string msg = "Exception caught: " + std::string(e.what());
		mexErrMsgTxt(msg.c_str());
	}

} // close mexFunction
