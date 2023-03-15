/**
 * @author Zorah Lähner (laehner@in.tum.de)
 * @version 03.2016
 */

#ifndef manifoldDijkstraLG_HPP
#define manifoldDijkstraLG_HPP

#define PRUNE false
#define DEBUG_MANIFOLD_DIJKSTRA false

#include <vector>
#include <list>
#include <chrono>
#include "MinHeapCG.hpp"
#include "arap.hpp"
#include "mex.h"
inline void updateArrays(      const long currentIdxInLayer,
                               const long currentLayer,
                               const long current3dSrc,
                               const long current3dTrgt,
                               const long nextIndex, 
                               const long nextIdxInLayer,
                               const long next2dSrc,
                               const long next2dTrgt,
                               const long next3dSrc,
                               const long next3dTrgt,
                               const double currentEnergy,
                               const double upperBound,
                               long* predecessors,
                               MinHeap* nextHeap,
                               double* result,
                               const double* all_energies,
                               const double* vertices,
                               const long nVertices,
                               const int dimVertices,
                               const double* triangles,
                               const long nTriangles,
                               const double* contour,
                               const long nContour,
                               const int dimContour,
                               const int nEdges3d,
                               const vector<std::tuple<long,long>>& edges3d,
                               const std::vector<Eigen::Quaterniond>& precRotations) {

    if (DEBUG_MANIFOLD_DIJKSTRA && nextIndex >= 6 * edges3d.size() * nContour) {
        std::cout << "result access beyond bounds " << nextIdxInLayer << " "<< edges3d.size() << std::endl;
        return;
    }
    if (result[nextIndex] < 0) {
        // index is not minimal yet
        double oldValue = nextHeap->peakKey(nextIdxInLayer);
        double newValue = 0;
        
        // handle linegraphPN overload (3d src might be just indicator how to enter 3d deg node)
        const long n3dSrc = nextIdxInLayer >= 4 * nEdges3d ? next3dTrgt : next3dSrc;
        const long c3dSrc = currentIdxInLayer >= 4 * nEdges3d ? current3dTrgt : current3dSrc;
        
        const long nProductNodesPerLayer = 6 * edges3d.size();
        newValue = currentEnergy + getArapCost(  vertices, nVertices, dimVertices, contour, nContour, dimContour,
                currentIdxInLayer, nextIdxInLayer, currentLayer, next2dSrc, next2dTrgt, current3dSrc,
                current3dTrgt, next3dSrc, next3dTrgt, c3dSrc, n3dSrc, nProductNodesPerLayer, precRotations);
        
        if (newValue <= upperBound && newValue < oldValue) {
            // new way is faster
            nextHeap->decrease(nextIdxInLayer, newValue);
            if (DEBUG_MANIFOLD_DIJKSTRA && nextIndex >= 6 * edges3d.size() * nContour) {
                std::cout << "predecessors access beyond bounds " << std::endl;
                return;
            }
            predecessors[nextIndex] = currentIdxInLayer + getLGStartLayerIndex(currentLayer, nVertices, nEdges3d);
        }
        
    }
}


void manifoldDijkstra(  const double* all_energies,
        const double upperbound,
        const LGProdNode* nodesPerLayer,
        const vector<list<long>* >&   mesh,
        const std::vector<std::list<std::tuple<long, long>>*>& mesh2, 
        const vector<std::tuple<long,long>>& edges3d, 
        const   double*                 vertices,
        const long                      nVertices,
        const int                       dimVertices,
        const   double*       triangles,
        const long            nTriangles,
        const   double*       contour,
        const long            nContour,
        const int             dimContour,
        double*               result,
        long*                 predecessors,
        const   int*          cut,
        const   int           cutID, // dijkstra will start only on vertices where cut[i]==cutID
        const std::vector<Eigen::Quaterniond>& precRotations) {
    std::cout << "  ManifoldDijkstra Iteration [";
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    using std::chrono::nanoseconds;
    auto t00 = high_resolution_clock::now();

    const long nEdges3d = edges3d.size();
    const long nEdges3dX2 = nEdges3d * 2;
    const long nEdges3dX4 = nEdges3d * 4;
            
    // initialise and find upper bound
    
    MinHeap** heapStruct = new MinHeap*[nContour];
    
    //const int numNodesOnLayer = 4 * nEdges3d + nVertices;
    const int numNodesOnLayer = 6 * nEdges3d;
    heapStruct[0] = new MinHeap(numNodesOnLayer);
    
    for (long i = numNodesOnLayer-1; i >= 0; i--) {
        if(cut[i] == cutID) {
        	heapStruct[0]->push(0, i);
        }
    }
    
    for (long i = 0; i < nContour; i++) {
        if (i > 0) {
            heapStruct[i] = new MinHeap(numNodesOnLayer);
        }
        
        for (long j=0; j < numNodesOnLayer; j++) {
            if (i > 0) {
                heapStruct[i]->push(std::numeric_limits<double>::infinity(), j);
            }
            
            result[i * numNodesOnLayer + j] = -1;
            predecessors[i * numNodesOnLayer + j] = -1;
        }
        
    }
    
    // run dijkstra
    long nextIndex;
    int level = 0;
    MinHeap* heap = heapStruct[level];
    if (PRUNE && DEBUG_MANIFOLD_DIJKSTRA) {
        std::cout << "LineGraph is pruned (no steeringAngles on 3D mesh with more than 90deg allowed)" << std::endl;
    }
    
    const float progressStepWidth = 1/20.0;
    float progressStep = 0;
    
    bool anyEnergyLowerThanUpperBound = false;
    while(!heap->isEmpty()) {
        if ( ((float)level/nContour) > progressStep) {
            mexPrintf("=");
            // workaround so progress bar actually is visible
            // https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
            mexEvalString("pause(.00001);");
            progressStep += progressStepWidth; 
        }

        pair<double, long> current = heap->pop();
        const int currentIdxInLayer = current.second;
        
        double currentEnergy = current.first;
        result[level * numNodesOnLayer + currentIdxInLayer] = currentEnergy;
        
        if(currentEnergy <= upperbound) {
            anyEnergyLowerThanUpperBound = true;
            //LGProdNode currentPN = getProductNodeInLayer(currentIdxInLayer, level, edges3d, nVertices);
            LGProdNode currentPN = nodesPerLayer[currentIdxInLayer];//getProductNodeInLayer(currentIdxInLayer, level, edges3d, nVertices);
            currentPN.idxSrc2d += level;
            currentPN.idxTrgt2d += level;
            
            /* 
             * current node == degenerate2d edge   => connect to non-deg this layer and deg2d this layer
             *              possible pruning: only connect to deg2d so that line    
             * current node == degenerate3d edge   => connect to non-deg next layer and deg3d next layer
             *              
             * current node == non-degenerate edge => connect to everything next layer
             *
             *
             */
            const bool isNonDeg = currentIdxInLayer < 2 * nEdges3d;
            const bool isDeg2d = (currentIdxInLayer < 4 * nEdges3d) && !isNonDeg;
            const bool isDeg3d = !isDeg2d && !isNonDeg;
            
            const long next2dSrc = currentPN.idxTrgt2d; 
            const long next3dSrc = currentPN.idxTrgt3d;
            // this is very confusing i know: there are nContour many layers but we need to acess until nContour+1
            // this is why we ensure the sufficient size of inputs in sm_2d_3d_dijkstra
            const bool isNotLastLayer = next2dSrc <= nContour;
            const bool isNotPreLastLayer = next2dSrc < nContour;
            const long nextLgStartLayerIdx = getLGStartLayerIndex(next2dSrc, nVertices, nEdges3d);
            
            // optimized:
            for (auto it : *mesh2.at(next3dSrc)) {
                const long next3dTrgt = std::get<0>(it);
                if (next3dTrgt == currentPN.idxSrc3d) {
                    // no going back i.e. we do not allow going current 3d target be vi if previous source was vi
                    continue;
                }
                const long edgeIdx = std::get<1>(it);
                const long nextLgStartLayerIndexForLoop = getLGStartLayerIndex(next2dSrc, nVertices, nEdges3d);
                
                // move to non-deg this or next layer 
                if ( isNotPreLastLayer && (isNonDeg || isDeg3d || isDeg2d)) {
                    const long next2dTrgt = next2dSrc + 1;

                    const long nextIdxInLayer = edgeIdx;
                    const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                    
                    updateArrays(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                              next2dSrc, next2dTrgt, next3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors,
                                              heapStruct[next2dSrc], result, all_energies, vertices, nVertices, dimVertices, triangles, 
                                              nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations);
                }
                
                // move to deg2d this layer
                if (isDeg2d) {
                    const long next2dTrgt = next2dSrc;
                     
                    const long nextIdxInLayer = getLGIndexInLayerForDeg2d(edgeIdx, nEdges3d);
                    const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                    updateArrays(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                              next2dSrc, next2dTrgt, next3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors,
                                              heapStruct[next2dSrc], result, all_energies, vertices, nVertices, dimVertices, triangles, 
                                              nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations);
                }
                
                // move to deg2d next layer
                if (isNotPreLastLayer && (isNonDeg)) {
                    const long next2dTrgt = next2dSrc;
                    
                    const long nextIdxInLayer = getLGIndexInLayerForDeg2d(edgeIdx, nEdges3d);
                    const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                    updateArrays(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                              next2dSrc, next2dTrgt, next3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors,
                                              heapStruct[next2dSrc], result, all_energies, vertices, nVertices, dimVertices, triangles, 
                                              nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations);
                }
            }
                
            // move to deg3d this or next layer
            if (isNotPreLastLayer && (isNonDeg || isDeg3d)) {
                const long next2dTrgt = next2dSrc + 1;
                const long next3dTrgt = next3dSrc;
                const long n3dSrc = currentPN.idxSrc3d;

                const long nextIdxInLayer = currentIdxInLayer >= 4 * nEdges3d ? 
                                                currentIdxInLayer : 
                                                currentIdxInLayer >= 2 * nEdges3d ? 
                                                    currentIdxInLayer + nEdges3dX2 : 
                                                    currentIdxInLayer + nEdges3dX4;
                const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                updateArrays(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                          next2dSrc, next2dTrgt, n3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors,
                                          heapStruct[next2dSrc], result, all_energies, vertices, nVertices, dimVertices, triangles, 
                                          nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations);                
            }
            
           
               
        }
        
        if(heap->isEmpty() && level < nContour - 1) {
            // early stopping since we cannot make progress anymore
            if (!anyEnergyLowerThanUpperBound) {
                if (DEBUG_MANIFOLD_DIJKSTRA) 
                    std::cout << "Early stopping bc no paths < upperbound found" << std::endl;

                // complete progress bar
                for (int l = level; l < nContour; l++) {
                    if ( ((float)level/nContour) > progressStep) {
                        std::cout << "=";
                        progressStep += progressStepWidth; 
                    }
                }
                break;
            }
            level++;
            heap = heapStruct[level];
        }
        
    }
    
    for(int i = 0; i < nContour; i++) {
        delete heapStruct[i];
    }
    


    auto t01 = high_resolution_clock::now();
    auto ms_____int = duration_cast<milliseconds>(t01 - t00);
    std::cout << "]  took: " << ms_____int.count() << "ms "<< std::endl;
    delete[] heapStruct;
}

#endif /* manifoldDijkstraLG_HPP */