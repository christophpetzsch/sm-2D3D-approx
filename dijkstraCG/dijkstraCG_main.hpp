/*
    Product Graph:
        @version 03.2016 author Zorah Lähner (laehner@in.tum.de)
    Adaption for Conjugate Graph:
        @version 03.2023 author Paul Rötzer (paul.roetzer@uni-bonn.de)
*/
#define DEBUG_DIJKSTRA_MAIN false 
#include "arap.hpp"

double dijkstra_main( 	const double* 	vertices,
                        const long 		nVertices,
                        const int 		dimVertices,
                        const double* 	triangles,
                        const long 		nTriangles,
                        const double* 	contour,
                        const long 		nContour,
                        const int 		dimContour,
                        std::vector<std::tuple<long,long>>& matching) {
    list<long> vPath;
    
    // build mesh data structure
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
    const long nProductNodes = nProductNodesPerLayer * nContour;
    const long firstIdxInLastLayer = nProductNodesPerLayer * (nContour-1);
    

    LGProdNode nodesPerLayer[nProductNodesPerLayer];
    for (int i = 0; i < nProductNodesPerLayer; i++) {
        nodesPerLayer[i] = getProductNodeInLayer(i, 0, edges3d, nVertices);
    }
	double* const result = new double[nProductNodes];
	long* const predecessors = new long[nProductNodes];

    std::vector<Eigen::Quaterniond> precRotations;
    if (PRECOMPUTE_ROTATIONS) {
        precRotations.reserve(nProductNodes);
        for (int i = 0; i < nContour; i++) {
            for (int j = 0; j < nProductNodesPerLayer; j++) {
                const long idx = nContour * nProductNodesPerLayer + j;
                LGProdNode cn = nodesPerLayer[j];
                // add level
                cn.idxSrc2d += i; cn.idxTrgt2d += i;
                long c2dTrgt = cn.idxTrgt2d;
                long c2dSrc = cn.idxSrc2d;
                c2dSrc =  cn.idxSrc2d == cn.idxTrgt2d  ? cn.idxSrc2d-1 : cn.idxSrc2d;
                if (c2dSrc < 0) c2dSrc = nContour-2;
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

    double* upperBoundEnergyArray = new double[nVertices];
	double* const all_energies = new double[nVertices * nContour];

	for (int i=0; i<nContour; ++i) {
		for (int j = 0; j < nVertices; ++j) {
			if (i == 0) {
                upperBoundEnergyArray[j] = 0.;
            }
            else {
                const int edgeIdx = std::get<1>(mesh2.at(j)->front()); 
                const int idx3dSrc = std::get<0>(mesh2.at(j)->front()); // we need this for LG overload
                const int idxInLayer = 4 * nEdges3d + edgeIdx;

                upperBoundEnergyArray[j] += getArapCost(vertices, nVertices, dimVertices, contour, nContour, dimContour, 
                            idxInLayer, idxInLayer, i-1, i, i+1, idx3dSrc, j, idx3dSrc, j, j, j, nProductNodesPerLayer,precRotations); 
                   
            }
        }
	}


    double upperbound = std::numeric_limits<double>::infinity();
    int upperBoundIdx = 0;
    for (int j = 0; j < nVertices; ++j) {
        if (upperBoundEnergyArray[j] < upperbound) {
            upperbound = upperBoundEnergyArray[j];
            upperBoundIdx = j;
        }
    }
    delete[] upperBoundEnergyArray;


	bool optimum = false;
	// first cut is the whole shape
	int* const cuts 		= new int[nVertices];
    int* cuts2 = new int[nProductNodesPerLayer];
	MinHeap upperBounds(nProductNodesPerLayer);

    for(int i=0; i < nProductNodesPerLayer; i++) {
        if (i < nVertices)
            cuts[i] = 0;
		cuts2[i] = 0;
	}
	int cutID = 0;
    int maxID = 0;
	upperBounds.push(0, 0);

	long meshIndex, minIndex;
	long pathBack, pathFront;
	double min;
    int noPathFoundCounter = 0; 
    std::vector<long> minIndices; minIndices.reserve(1000);
    std::set<int> finishedCuts;
    std::set<int> upperboundCuts;
    std::vector<long> upperboundPath;

	// keep searching until there are no more cuts or the optimum (a loop) was found
    int loopCounter = 0;
    int iterOptFound = -1;
    double gapOptFound = upperbound - upperBounds.peak();
	while(!upperBounds.isEmpty() && !optimum) {
        std::cout << "Upperbound: " << upperbound << " Lowerbound: " << 
                        upperBounds.peak() << " gap = "<< upperbound - upperBounds.peak() << std::endl; 
        loopCounter++;
		cutID = upperBounds.pop().second;

		// run dijkstra on product manifold
		manifoldDijkstra(all_energies, upperbound, nodesPerLayer, mesh, mesh2, edges3d, vertices, nVertices, 
                         dimVertices, triangles, nTriangles, contour, nContour, 
                         dimContour, result, predecessors, cuts2, cutID, precRotations);

		// find minimum
		min = std::numeric_limits<double>::infinity();//result[firstIdxInLastLayer];
		minIndex = -1;//firstIdxInLastLayer;
        minIndices.clear();
		for (int i = firstIdxInLastLayer; i < nProductNodes; i++) {
			const long idxInLayer = i - firstIdxInLastLayer;
			if (cuts2[idxInLayer] == cutID && result[i] < min && result[i] != -1) {
				min = result[i];
				minIndex = i;
                minIndices.push_back(minIndex - firstIdxInLastLayer);
			}
            if (cuts2[idxInLayer] == cutID && result[i] == -1) {
                // rule out vertices which definitely cannot be on an optimal path anymore
                cuts2[idxInLayer] = -1;
            }
		}
        if (DEBUG_DIJKSTRA_MAIN) std::cout << "#minIndices " << minIndices.size() << std::endl;
        if (minIndex == -1) {
            finishedCuts.insert(cutID);
            // no path for this shape cut found
            std::cout << "  No path for this shapeCut found. Looking into another. Shape cuts left: " << upperBounds.size() << std::endl;
            noPathFoundCounter++;
            if (noPathFoundCounter > 100) // make sure we do not spin forever
                break;
            continue;
        }
        noPathFoundCounter = 0; // reset to avoid non-desired returns

		// create path
		vPath.clear();
        calcPath(minIndex, predecessors, vPath);

        
        pathFront = vPath.front();
        pathBack = vPath.back()  - firstIdxInLastLayer;
        if (DEBUG_DIJKSTRA_MAIN) std::cout << "pathFront:" << pathFront  << " pathBack:" << pathBack << std::endl;
		if(!vPath.empty() && pathFront == pathBack) {
            if (DEBUG_DIJKSTRA_MAIN && min > upperBounds.peak()) std::cout <<  "min > upperBounds.peak()" << std::endl;
			if (min <= upperBounds.peak()) {
                if (DEBUG_DIJKSTRA_MAIN) std::cout <<  "Optimum found" << std::endl;
				optimum = true;
                if (DEBUG_DIJKSTRA_MAIN) {
                    if (min == upperbound) {
                        std::cout << "Found optimum in last iteration" << std::endl;
                    }
                    else {
                        std::cout << "Found optimum " << loopCounter-iterOptFound << 
                                " iters ago but couldnt identify as such." << std::endl;
                        std::cout << "Gap when optimum found " << gapOptFound << std::endl;
                    }
                
                }
			} 
            else {
                upperboundCuts.insert(cutID);
				upperBounds.push(min, cutID);
                if (min <= upperbound) {
                    upperboundPath.assign(vPath.begin(), vPath.end());
                    upperbound = min;
                    iterOptFound = loopCounter;
                    gapOptFound = upperbound - upperBounds.peak();
                }
			}
		}  
        else {
			// cut shape in two smaller parts
			maxID++;
            
            for (auto it: minIndices) {
                const int minIdxOnLastLayer = it + firstIdxInLastLayer;
                vPath.clear();
                calcPath(minIdxOnLastLayer, predecessors, vPath);
                const long currentPathFront = vPath.front();
                
                if (it == currentPathFront) {
                    if (result[minIdxOnLastLayer] <= upperbound) {
                        upperbound = result[it + firstIdxInLastLayer];
                        upperboundPath.assign(vPath.begin(), vPath.end());
                        iterOptFound = loopCounter;
                        gapOptFound = upperbound - upperBounds.peak();
                        if (DEBUG_DIJKSTRA_MAIN) std::cout << "Found new upperbound in min indices " << upperbound  << std::endl;
                    }
                }
            }
            shapeCut2(edges3d, mesh2, nodesPerLayer, vertices, nVertices, 
                      dimVertices, nProductNodesPerLayer, 
                      cuts, cuts2, maxID, pathFront, pathBack);


            if (DEBUG_SHAPE_CUT ) {
                int shapeCutsCount[maxID+1];
                for (int j = 0; j< maxID+1; j++ ) {
                    shapeCutsCount[j] = 0;
                }
                for (int i = 0; i < nProductNodesPerLayer; i++) {
                    if (cuts2[i] <= maxID && cuts2[i] >= 0)
                        shapeCutsCount[cuts2[i]] += 1;
                }
                std::cout << "Active cuts: " << std::endl;
                for (int j = 0; j< maxID+1; j++ ) {
                    auto it1 = finishedCuts.find(j); const bool isFinishedCut = it1 != finishedCuts.end();
                    if (isFinishedCut) continue;
                    auto it2 = upperboundCuts.find(j); const bool isUpperbound = it2 != upperboundCuts.end();
                    std::cout << "Cut " << j << " " << shapeCutsCount[j]<< " nodes " ;
                    //if (isFinishedCut) std::cout << "- finished";
                    if (isUpperbound) std::cout << "- valid path found (upperbound cut)";
                    std::cout << std::endl;
                }
            }

			upperBounds.push(min, cutID);
			upperBounds.push(min, maxID);
		}
	}

    const double energy = minIndex == -1 ? upperbound : result[minIndex];
    // if minIndex==-1 set minIndex to upperbound path 
    if (minIndex == -1) {
        if (upperboundPath.size() == 0) {
            const long idx3d = upperBoundIdx;
            for (long j = 0; j < nContour; j++)
                matching.push_back(std::tuple<long, long>{idx3d, j});
            return energy;
        }
        // return best suboptimal path
        vPath.assign(upperboundPath.begin(), upperboundPath.end());
    }

    // clean memory
	for (int i=0; i<nVertices; i++) {
		delete mesh.at(i);
        delete mesh2.at(i);
	}

	delete[] result;
	delete[] predecessors;
	delete[] cuts;
    delete[] cuts2;
	delete[] all_energies;
    
    
    // compute matching from path
    for (auto const& nodeId: vPath) {
        LGProdNode pn = getProductNode(nodeId, edges3d, nVertices);
        const long idx3d = getLGSrc3dVertex(nodeId, edges3d, nVertices);
        const long idx2d = pn.idxSrc2d % (nContour-1);
        
        const double dist2Other2d1 = getLocalThickness2d(idx2d, contour, dimContour);
        const double dist2Other2d2 = getLocalThickness2d(pn.idxTrgt2d, contour, dimContour);
        const double dist2Other3d1 = getLocalThickness3d(idx3d, vertices, dimVertices);
        const double dist2Other3d2 = getLocalThickness3d(pn.idxTrgt3d, vertices, dimVertices);
        if (DEBUG_DIJKSTRA_MAIN) std::cout << idx2d << "-" << idx3d  << 
                " >>  " << pn.idxTrgt2d % (nContour-1) << "-" << pn.idxTrgt3d << " " << std::endl;
        
        matching.push_back(std::tuple<long, long>{idx3d, idx2d});
    }

	return energy;
}
