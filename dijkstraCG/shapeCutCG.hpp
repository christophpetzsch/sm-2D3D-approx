/**
		@author Zorah Lähner (laehner@in.tum.de)
    @version 03.2016
*/

#include "MinHeapCG.hpp"
#define DEBUG_SHAPE_CUT false

using std::vector;
using std::list;

// dijkstra from two points cutting the shape in half
// more or less the voronoi cells around the points
void shapeCut(	const 	vector<list<long>* >& 	mesh,
				 	long 					nVertices,
					double* 			result,
					int* 					cut,
				 	int 					newId,
				 	long 					vertex0,
				 	long 					vertex1
			 ) {

	const int oldId = cut[vertex0];
	if (oldId != cut[vertex1]) {
		return;
    }


    // initialise heap
    MinHeap heap(nVertices);
	heap.push(0, vertex1);
	heap.push(0, vertex0);

	for (long i = 0; i < nVertices; i++) {
		if(cut[i] != oldId) {
			result[i] = 0;
		} 
        else {
			result[i] = -1;
            heap.push(std::numeric_limits<double>::infinity(), i);
		}
	}

	cut[vertex1] = newId; // all verticies close to vertex1 will end up with newId

	while(!heap.isEmpty()) {
		pair<double, long> current = heap.pop();
		result[current.second] = current.first;

		double oldValue, newValue;

		for(auto it = mesh.at(current.second)->begin(); it != mesh.at(current.second)->end(); it++) {
			if(result[*it] < 0) {
				// index is not minimal yet
				oldValue = heap.peakKey(*it);
				newValue = current.first + 1;
				if (newValue < oldValue) {
					// new way is faster
					heap.decrease(*it, newValue);
					cut[*it] = cut[current.second];
				}
			}
		}
	}

}


int getEdgeIdxFromLGIndex(const int lgIndex, const int nEdges3d) {
    int edgeIdx = lgIndex;
    if (lgIndex >= 4 * nEdges3d)
        edgeIdx =  lgIndex - 4 * nEdges3d;
    else if (lgIndex >= 2 * nEdges3d)
        edgeIdx =  lgIndex - 2 * nEdges3d;
    return edgeIdx;
}

int getEdgeIdxFromEIdx(const int edgeIdx, const int nEdges3d) {
    if (edgeIdx > nEdges3d)
        return edgeIdx - nEdges3d;
    return edgeIdx;
}

int getFirstVertexIdx(const int lgIndex, const std::vector<std::tuple<long,long>>& edges3d) {
    const int nEdges3d = edges3d.size();
    int edgeIdx = lgIndex;
    if (lgIndex >= 4 * nEdges3d)
        edgeIdx =  lgIndex - 4 * nEdges3d;
    else if (lgIndex >= 2 * nEdges3d)
        edgeIdx =  lgIndex - 2 * nEdges3d;
    
    if (edgeIdx >= nEdges3d)
        return std::get<1>(edges3d.at(edgeIdx-nEdges3d));
    else
        return std::get<0>(edges3d.at(edgeIdx));
}

int getSecondVertexIdx(const int lgIndex, const std::vector<std::tuple<long,long>>& edges3d) {
    const int nEdges3d = edges3d.size();
    int edgeIdx = lgIndex;
    if (lgIndex >= 4 * nEdges3d)
        edgeIdx =  lgIndex - 4 * nEdges3d;
    else if (lgIndex >= 2 * nEdges3d)
        edgeIdx =  lgIndex - 2 * nEdges3d;
    
    if (edgeIdx >= nEdges3d)
        return std::get<0>(edges3d.at(edgeIdx-nEdges3d));
    else
        return std::get<1>(edges3d.at(edgeIdx));
}


void shapeCut2(	const std::vector<std::tuple<long,long>>& edges3d,
                const std::vector<std::list<std::tuple<long, long>>*>& mesh2, 
                const LGProdNode* nodesPerLayer,
                const double*   features3d,
                const int       nVertices3d,
                const int       dimFeatures3d,
                long 			numNodesPerLayer,
                int* 			cut,
                int* 			cut2,
                int 			newId,
                long 			nodeId0,
                long 			nodeId1) {
    
    if (nodeId0 == nodeId1 && DEBUG_SHAPE_CUT) {
        std::cout << "Optimal Values Found... what you are searching here??" << std::endl;
        return;
    }
	const int oldId = cut2[nodeId0];
    if (DEBUG_SHAPE_CUT) std::cout << " Old ID " << oldId << std::endl;
	if (oldId != cut2[nodeId1] && DEBUG_SHAPE_CUT) {
        std::cout << "sorrry, but you should not input this here..." << std::endl;
		return;
    }
    
    cut2[nodeId1] = newId;
    
    const int nEdges3d = edges3d.size();
    double edgeLengths[nEdges3d];
    for (int i = 0; i < nEdges3d; i++) {
        double edgeLength = 0; 
        double edge[3];
        const int src = std::get<0>(edges3d.at(i));
        const int trgt = std::get<1>(edges3d.at(i));
        for (int j = 0; j < 3; j++) {
            const double srcD = features3d[dimFeatures3d * src + j];
            const double trgtD = features3d[dimFeatures3d * trgt + j];
            edge[j] = trgtD - srcD;
            edgeLength += edge[j] * edge[j];
        }
        edgeLengths[i] = std::sqrt(edgeLength);
    }
    
    
    // 0) find start vertices
    int vertex0 = getSecondVertexIdx(nodeId0, edges3d);
    int vertex1 = getSecondVertexIdx(nodeId1, edges3d);
    
    /*
     * vertex1 != vertex0
     *
     */
    if (vertex1 != vertex0) {
        // 1) Cut shape in half as before
        double result[nVertices3d];
        MinHeap heap(nVertices3d);
        heap.push(0, vertex1);
        heap.push(0, vertex0);

        for (long i = 0; i < nVertices3d; i++) {
            if (cut[i] != oldId)
                result[i] = 0;
            else {
                result[i] = -1;
                heap.push(std::numeric_limits<double>::infinity(), i);
            }
        }

        cut[vertex1] = newId; // all verticies close to vertex1 will end up with newId

        while(!heap.isEmpty()) {
            pair<double, long> current = heap.pop();
            result[current.second] = current.first;

            double oldValue, newValue;
            for (auto it : *mesh2.at(current.second)) {
                const long nextVertex = std::get<0>(it);
                const long edgeIdx = getEdgeIdxFromEIdx(std::get<1>(it), nEdges3d);
                if (result[nextVertex] < 0) {
                    oldValue = heap.peakKey(nextVertex);
                    newValue = current.first + edgeLengths[edgeIdx];
                    if (newValue < oldValue) {
                        heap.decrease(nextVertex, newValue);
                        cut[nextVertex] = cut[current.second];
                    }
                }
            }
        }

        // 2) transfer cut to product graph
        int nodesOnNewCut = 0;
        for (long i = 0; i < numNodesPerLayer; i++) {
            if (cut2[i] == oldId) {
                const LGProdNode pn = nodesPerLayer[i];
                const long src3d = pn.idxSrc3d;
                const long trgt3d = pn.idxTrgt3d;

                if (cut[trgt3d] == newId) {
                    cut2[i] = newId;
                    nodesOnNewCut++;
                }
            }
        }

        if (nodesOnNewCut == 0 && DEBUG_SHAPE_CUT) {
            std::cout << "No nodes on new cut: vertex0 " << vertex0 << " " << "vertex1 " << vertex1 <<std::endl;
        }
    
    
    }
    
    /*
     * trgt1 == trgt0
     *
     */
    else {
        int src0 = getFirstVertexIdx(nodeId0, edges3d);
        int src1 = getFirstVertexIdx(nodeId1, edges3d);
         /*
         * src0 == src1
         *
         */
        if (src0 == src1) {
            // just move this one vertex to new cut and the other one containing
            // same source and target and are not same as first idx
            long otherId = -1;
            if (nodeId0 < 2 * nEdges3d) {
                if (nodeId1 < 2 * nEdges3d) {
                    if (DEBUG_SHAPE_CUT) std::cout << "Cannot happen nodeId0 < 2 * nEdges3d && nodeId1 < 2 * nEdges3d" << std::endl; 
                }
                else if (nodeId1 < 4 * nEdges3d) {
                    otherId = nodeId0 + 4 * nEdges3d;
                }
                else {
                    otherId = nodeId0 + 2 * nEdges3d;
                }
            }
            else if (nodeId0 < 4 * nEdges3d) {
                if (nodeId1 < 2 * nEdges3d) {
                    otherId = nodeId0 + 4 * nEdges3d;
                }
                else if (nodeId1 < 4 * nEdges3d) {
                    if (DEBUG_SHAPE_CUT) std::cout << "Cannot happen nodeId0 < 4 * nEdges3d && nodeId1 < 4 * nEdges3d" << std::endl; 
                }
                else {
                    otherId = nodeId0 - 2 * nEdges3d;
                }
            }
            else {
                if (nodeId1 < 2 * nEdges3d) {
                    otherId = nodeId0 - 2 * nEdges3d;
                }
                else if (nodeId1 < 4 * nEdges3d) {
                    otherId = nodeId0 - 4 * nEdges3d;
                }
                else {
                    if (DEBUG_SHAPE_CUT) std::cout << "Cannot happen nodeId0 >= 4 * nEdges3d && nodeId1 >= 4 * nEdges3d" << std::endl; 
                }
            }
            if (otherId == -1 && DEBUG_SHAPE_CUT)
                std::cout << "Cannot happen, didnt find other id" << std::endl; 
            if (cut2[otherId] == oldId)
                cut2[otherId] = newId;
            return;
        }
        /*
         * src0 != src1
         *
         */
        else {
            if (DEBUG_SHAPE_CUT)
                std::cout << "Entering 3d vertex from opp direction than than nodeId0" << std::endl; 
            // move all product nodes related to src1 ---> trgt1 edge of 3d shape to new cut
            int otherIds[2];
            if (nodeId1 < 2 * nEdges3d) {
                otherIds[0] = nodeId1 + 2 * nEdges3d;
                otherIds[1] = nodeId1 + 4 * nEdges3d;
            }
            else if (nodeId1 < 4 * nEdges3d) { 
                otherIds[0] = nodeId1 - 2 * nEdges3d;
                otherIds[1] = nodeId1 + 2 * nEdges3d;
            }
            else {
                otherIds[0] = nodeId1 - 4 * nEdges3d;
                otherIds[1] = nodeId1 - 2 * nEdges3d;
            }
            for (int i = 0; i < 2; i++) {
                if (cut2[otherIds[i]] == oldId)
                    cut2[otherIds[i]] = newId;
            }
        }
    }
    return;
    
    
    
    
    
    
    /*if (vertex1 == vertex0) {
        std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Same Target"<< std::endl;
        if (nodeId0 >= 2 * nEdges3d) {
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Current pair contains deg2d";
        }
        else if (nodeId0 >= 4 * nEdges3d){
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Current pair contains deg3d";
        }
        else {
            std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  Current pair contains nondeg";
        }
        if (nodeId1 >= 2 * nEdges3d){
            std::cout << " and deg2d"<< std::endl;
        }
        else if (nodeId1 >= 4 * nEdges3d){
            std::cout << " and deg3d"<< std::endl;
        }
        else {
            std::cout << " and nondeg"<< std::endl;
        }*/
        
        
        
        
        /*if (DEBUG_SHAPE_CUT) std::cout << "Both vertices identical, using fallback" << std::endl;
        // move all product nodes which are connected to vertex0 to other cut except nodeId0
        for (auto it : *mesh2.at(vertex0)) {
            const long nextVertex = std::get<0>(it);
            long edgeIdx = getEdgeIdxFromEIdx(std::get<1>(it), nEdges3d);
            // only move incoming so turn around edge (by using opposite orientation of edgeIdx)
            if (edgeIdx >= nEdges3d) 
                edgeIdx = edgeIdx - nEdges3d;
            else
                edgeIdx = edgeIdx - nEdges3d;
            
            
			for (int i = 0; i < 3; i+=2 ) {
                const int nodeID = edgeIdx + i * nEdges3d;
                if (edgeIdx == nodeId0)
                    continue;
                if (cut2[nodeID] == oldId) {
                    cut2[nodeID] = newId;
                }
			}
		}
        return;*/
        
        
        /*if (DEBUG_SHAPE_CUT) std::cout << "Both vertices identical, switching" << std::endl;
        vertex1 = getFirstVertexIdx(nodeId1, edges3d);
        if (vertex0 == vertex1) {
            if (DEBUG_SHAPE_CUT) std::cout << "Switching again" << std::endl;
            vertex0 = getFirstVertexIdx(nodeId0, edges3d);
            vertex1 = getSecondVertexIdx(nodeId1, edges3d);
        }
    }
    
    if (cut[vertex1] != oldId) {
        if (DEBUG_SHAPE_CUT) std::cout << "cut[vertex1] != oldID switching" << std::endl;
        vertex1 = getFirstVertexIdx(nodeId1, edges3d);
        if (cut[vertex1] != oldId) {
            if (DEBUG_SHAPE_CUT) std::cout << "Didnt work... fallback: splitting cut in two equal parts" << std::endl;
            std::vector<int> nodesInCut; nodesInCut.reserve(numNodesPerLayer/10);
            for (int i = 0; i < numNodesPerLayer; i++) {
                if (cut2[i] == oldId) {
                    nodesInCut.push_back(i);
                }
            }
            if (DEBUG_SHAPE_CUT && nodesInCut.size() == 1)
                std::cout << "Only one node in cut, this must be wrong..." << std::endl;
            
            
            // split shape into two equal parts but make ensure that deg 2d and deg 3d of the same shape
            // are in the same part as long as possible
            const int offsets[] = {0, 2 * nEdges3d, 4*nEdges3d};
            const int numNodesInCurrentCut = nodesInCut.size();
            if (numNodesInCurrentCut == 1) {
                return; // we already have moved one vertex to the other cut
            }
            if (numNodesInCurrentCut == 2) {
                if (nodesInCut.at(0) != nodeId0 && cut2[nodesInCut.at(0)] == oldId)
                    cut2[nodesInCut.at(0)] = newId;
                else if (nodesInCut.at(1) != nodeId0 && cut2[nodesInCut.at(1)] == oldId)
                    cut2[nodesInCut.at(1)] = newId;
                else {
                    if (DEBUG_SHAPE_CUT) 
                        std::cout <<  "Cannot assign any vertex to new cut, this should not happen" << std::endl;
                }
                return;
            }
            int numNodesInNewCut = 1; // we already have nodeId1 in cut 
            const int stoppingNumber = std::floor(numNodesInCurrentCut/2);
            for (int i = 0; i < numNodesInCurrentCut; i++) {
                for (int j = 0; j < 3; j++) {
                    const int nodeId = nodesInCut.at(i) + offsets[j];
                    if (cut2[nodeId] == oldId) {
                        cut2[nodeId] = newId;
                        numNodesInNewCut++;
                    }
                    if ( numNodesInNewCut >= stoppingNumber)
                        break;
                }
                if (numNodesInNewCut >= stoppingNumber)
                    break;
            }
            return;
        }
    }
    
    if (DEBUG_SHAPE_CUT) std::cout <<  "vertex0 " << vertex0 << " vertex1 " << vertex1 << std::endl;
    
    if (vertex1 == vertex0 && DEBUG_SHAPE_CUT) {
        std::cout << "Both vertices identical, this will not work" << std::endl; 
        return;
    }
    
    // 1) Cut shape in half as before
    double result[nVertices3d];
    MinHeap heap(nVertices3d);
	heap.push(0, vertex1);
	heap.push(0, vertex0);

	for (long i = 0; i < nVertices3d; i++) {
		if (cut[i] != oldId)
			result[i] = 0;
        else {
			result[i] = -1;
            heap.push(std::numeric_limits<double>::infinity(), i);
		}
	}
    
	cut[vertex1] = newId; // all verticies close to vertex1 will end up with newId

	while(!heap.isEmpty()) {
		pair<double, long> current = heap.pop();
		result[current.second] = current.first;

		double oldValue, newValue;
        for (auto it : *mesh2.at(current.second)) {
            const long nextVertex = std::get<0>(it);
            const long edgeIdx = getEdgeIdxFromEIdx(std::get<1>(it), nEdges3d);
			if (result[nextVertex] < 0) {
				oldValue = heap.peakKey(nextVertex);
				newValue = current.first + edgeLengths[edgeIdx];
				if (newValue < oldValue) {
					heap.decrease(nextVertex, newValue);
					cut[nextVertex] = cut[current.second];
				}
			}
		}
	}
    
    // 2) transfer cut to product graph
    int nodesOnNewCut = 0;
    for (long i = 0; i < numNodesPerLayer; i++) {
        if (cut2[i] == oldId) {
            const LGProdNode pn = nodesPerLayer[i];
            const long src3d = pn.idxSrc3d;
            const long trgt3d = pn.idxTrgt3d;
            
            if (cut[trgt3d] == newId) {
                cut2[i] = cut2[i] = newId; // cut[src3d]; // 
                nodesOnNewCut++;
            }
        }
    }
    
    if (nodesOnNewCut == 0 && DEBUG_SHAPE_CUT) {
        std::cout << "No nodes on new cut: vertex0 " << vertex0 << " " << "vertex1 " << vertex1 <<std::endl;
    }*/

}
