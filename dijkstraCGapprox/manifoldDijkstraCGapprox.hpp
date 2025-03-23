/**
 * @author Zorah Lähner (laehner@in.tum.de)
 * @version 03.2016
 */

#ifndef manifoldDijkstraLG_HPP
#define manifoldDijkstraLG_HPP

#define PRUNE false
#define DEBUG_MANIFOLD_DIJKSTRA true

#include <vector>
#include <list>
#include <chrono>
#include "MinHeapCGapprox.hpp"
#include "arapapprox.hpp"
#include "mex.h"
inline void updateArrays_CG_CG(const long currentIdxInLayer,
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
                               long *predecessors_CG,
                               const long numNodesOnLayer,
                               MinHeap *nextHeap,
                               double *result_CG,
                               const double *all_energies,
                               const double *vertices,
                               const long nVertices,
                               const int dimVertices,
                               const double *triangles,
                               const long nTriangles,
                               const double *contour,
                               const long nContour,
                               const int dimContour,
                               const long nEdges3d,
                               const vector<std::tuple<long, long>> &edges3d,
                               const std::vector<Eigen::Quaterniond> &precRotations,
                               const int* landmarks_2d,
                               const int* landmarks_3d,
                               const int nLandmarks)
{
    if(next3dTrgt == current3dSrc) return;
    if (DEBUG_MANIFOLD_DIJKSTRA && nextIndex >= 6 * edges3d.size() * (nContour + 1))
    {
        std::cout << "result access beyond bounds " << nextIdxInLayer << " " << nextIndex - 6 * edges3d.size() * (nContour + 1) << " " << edges3d.size() << std::endl;
        return;
    }
    if (result_CG[getIndexForArrays(nextIndex, numNodesOnLayer)] < 0)
    {
        // index is not minimal yet
        double oldValue = nextHeap->peakKey(nextIdxInLayer);
        double newValue = 0;

        newValue = currentEnergy + getArapCost(vertices, nVertices, dimVertices, contour, nContour, dimContour,
                                               currentIdxInLayer, nextIdxInLayer, currentLayer, next2dSrc, next2dTrgt, current3dSrc,
                                               current3dTrgt, next3dSrc, next3dTrgt, numNodesOnLayer, precRotations);
        if(nLandmarks != 0)
        {
            if(landmarks_2d[1]==1)
            {
                if(landmarks_3d[0] != -1 && currentLayer == 0 && (next2dSrc!=1 || current3dSrc!=landmarks_3d[0] || current3dTrgt!=landmarks_3d[1])) newValue += 1000000;
            }
            else
            {
                if(landmarks_3d[0] != -1 && currentLayer == 0 && (currentIdxInLayer<4*nEdges3d || current3dTrgt != landmarks_3d[0])) newValue += 1000000;
            }
        }
        if (newValue <= upperBound && newValue < oldValue)
        {
            // new way is faster
            nextHeap->decrease(nextIdxInLayer, newValue);
            if (DEBUG_MANIFOLD_DIJKSTRA && nextIndex >= 6 * edges3d.size() * (nContour + 1))
            {
                std::cout << "predecessors access beyond bounds " << std::endl;
                return;
            }
            predecessors_CG[getIndexForArrays(nextIndex, numNodesOnLayer)] = currentIdxInLayer;
        }
    }
}

// Going from the CG node i1->i2 into the SG node i2, thus no new costs are calculated
inline void updateArrays_CG_SG(const long current3d,
                               const long edge,
                               const double currentEnergy,
                               const double upperBound,
                               long *entering_edge_SG,
                               const long numNodesOnLayer,
                               MinHeap *nextHeap,
                               double *result_SG,
                               const double *all_energies,
                               const double *vertices,
                               const long nVertices,
                               const int dimVertices,
                               const double *triangles,
                               const long nTriangles,
                               const double *contour,
                               const long nContour,
                               const int dimContour,
                               const long nEdges3d,
                               const vector<std::tuple<long, long>> &edges3d,
                               const std::vector<Eigen::Quaterniond> &precRotations,
                               const int* landmarks_2d,
                               const int* landmarks_3d,
                               const int nLandmarks)
{
    if (result_SG[current3d] < 0)
    {
        // index is not minimal yet
        double oldValue = nextHeap->peakKey(current3d);
        double newValue = currentEnergy;

        long reduced_edge = edge % (2 * nEdges3d);
        long local_vertex;
        if(reduced_edge >= nEdges3d)
            {
                local_vertex = std::get<0>(edges3d[reduced_edge-nEdges3d]);
            } else 
            {
                local_vertex = std::get<1>(edges3d[reduced_edge]);
        }
        if(local_vertex != current3d){
            std::cout << "wrong CG->SG code, " << local_vertex << " vs. " << current3d << std::endl;
        }

        if (newValue <= upperBound && newValue < oldValue)
        {
            // new way is faster
            nextHeap->decrease(current3d, newValue);
            entering_edge_SG[current3d] = edge;
        }
    }
}

// Going from the SG node i0 into the CG node i0->i1, thus costs are calculated
inline void updateArrays_SG_CG(const long first_edge_Src,
                               const long first_edge_Trgt,
                               const long second_edge_Src,
                               const long second_edge_Trgt,
                               const long previous2d,
                               const long current2d,
                               const long next2d,
                               const long first_edge,  // local CG index of edge going into (current2d, current3d)
                               const long second_edge, // local CG index of edge going from (current2d, current3d) to (next2d, next3d)
                               const double currentEnergy,
                               const double upperBound,
                               long *predecessors_CG,
                               const long numNodesOnLayer,
                               MinHeap *nextHeap,
                               double *result_CG,
                               const double *all_energies,
                               const double *vertices,
                               const long nVertices,
                               const int dimVertices,
                               const double *triangles,
                               const long nTriangles,
                               const double *contour,
                               const long nContour,
                               const int dimContour,
                               const long nEdges3d,
                               const vector<std::tuple<long, long>> &edges3d,
                               const std::vector<Eigen::Quaterniond> &precRotations,
                               const int* landmarks_2d,
                               const int* landmarks_3d,
                               const int nLandmarks)
{
    if(second_edge_Trgt == first_edge_Src) return;
    if (result_CG[getIndexForArrays(second_edge+current2d*numNodesOnLayer, numNodesOnLayer)] < 0)
    {
        // index is not minimal yet
        double oldValue = nextHeap->peakKey(second_edge);
        double newValue = 0;

        newValue = currentEnergy + getArapCost(vertices, nVertices, dimVertices, contour, nContour, dimContour, first_edge, second_edge, previous2d, current2d, next2d,
                                               first_edge_Src, first_edge_Trgt, second_edge_Src, second_edge_Trgt, numNodesOnLayer, precRotations);
        if(nLandmarks != 0)
        {
            if(landmarks_2d[1] == 1)
            {
                if(landmarks_3d[0] != -1 && second_edge_Trgt != landmarks_3d[0]) newValue += 1000000;
            }
            else
            {
                if(landmarks_3d[0] != -1 && second_edge_Trgt != landmarks_3d[0] || second_edge < 4*nEdges3d) newValue += 1000000;
            }
        }
        if (newValue <= upperBound && newValue < oldValue)
        {
            // new way is faster
            nextHeap->decrease(second_edge, newValue);
            predecessors_CG[getIndexForArrays(second_edge+current2d*numNodesOnLayer, numNodesOnLayer)] = first_edge;
        }
    }
}

inline void updateArrays_SG_SG(const long first_edge_Src,
                               const long first_edge_Trgt,
                               const long second_edge_Src,
                               const long second_edge_Trgt,
                               const long previous2d,
                               const long current2d,
                               const long next2d,
                               const long first_edge,  // local CG index of edge going into (current2d, current3d)
                               const long second_edge, // local CG index of edge going from (current2d, current3d) to (next2d, next3d)
                               const double currentEnergy,
                               const double upperBound,
                               long *entering_edge_SG,
                               const long numNodesOnLayer,
                               MinHeap *nextHeap,
                               double *result_SG,
                               const double *all_energies,
                               const double *vertices,
                               const long nVertices,
                               const int dimVertices,
                               const double *triangles,
                               const long nTriangles,
                               const double *contour,
                               const long nContour,
                               const int dimContour,
                               const long nEdges3d,
                               const vector<std::tuple<long, long>> &edges3d,
                               const std::vector<Eigen::Quaterniond> &precRotations,
                               const int* landmarks_2d,
                               const int* landmarks_3d,
                               const int nLandmarks)
{
    if(second_edge_Trgt == first_edge_Src) return;
    if (result_SG[second_edge_Trgt + (next2d - 2) * nVertices] < 0)
    {
        // index is not minimal yet
        double oldValue = nextHeap->peakKey(second_edge_Trgt);
        double newValue = 0;

        
        newValue = currentEnergy + getArapCost(vertices, nVertices, dimVertices, contour, nContour, dimContour, first_edge, second_edge, previous2d, current2d, next2d,
                                               first_edge_Src, first_edge_Trgt, second_edge_Src, second_edge_Trgt, numNodesOnLayer, precRotations);
        for(int i=1; i<nLandmarks; ++i)
        {
            if(landmarks_2d[i] == next2d && landmarks_3d[i] != second_edge_Trgt) newValue += 1000000;
        }
        if (newValue <= upperBound && newValue < oldValue)
        {
            // new way is faster
            nextHeap->decrease(second_edge_Trgt, newValue);
            entering_edge_SG[second_edge_Trgt + (next2d - 2) * nVertices] = second_edge;
        }
    }
}

void manifoldDijkstra(const double *all_energies,
                      const double upperbound,
                      const LGProdNode *nodesPerLayer,
                      const vector<list<long> *> &mesh,
                      const std::vector<std::list<std::tuple<long, long>> *> &mesh2,
                      const vector<std::tuple<long, long>> &edges3d,
                      const double *vertices,
                      const long nVertices,
                      const int dimVertices,
                      const double *triangles,
                      const long nTriangles,
                      const double *contour,
                      const long nContour,
                      const int dimContour,
                      double *result_CG,
                      double *result_SG,
                      long *predecessors_CG,
                      long *entering_edge_SG,
                      const int *cut,
                      const int cutID, // dijkstra will start only on vertices where cut[i]==cutID
                      const std::vector<Eigen::Quaterniond> &precRotations,
                      const int* landmarks_2d,
                      const int* landmarks_3d,
                      const int nLandmarks)
{
    std::cout << "  ManifoldDijkstra Iteration [";
    using std::chrono::duration;
    using std::chrono::duration_cast;
    using std::chrono::high_resolution_clock;
    using std::chrono::milliseconds;
    using std::chrono::nanoseconds;
    auto t00 = high_resolution_clock::now();

    const long nEdges3d = edges3d.size();
    const long nEdges3dX2 = nEdges3d * 2;
    const long nEdges3dX4 = nEdges3d * 4;

    // initialise and find upper bound

    MinHeap **heapStruct = new MinHeap *[nContour+1];

    const long numNodesOnLayer = 6 * nEdges3d;
    heapStruct[0] = new MinHeap(numNodesOnLayer);

    for (long i = numNodesOnLayer - 1; i >= 0; i--)
    {
        if (cut[i] == cutID)
        {
            heapStruct[0]->push(0, i);
        }
    }

    for (long i = 0; i < nContour + 1; i++)
    {
        if (i > 1 && i < nContour)
        {
            heapStruct[i] = new MinHeap(nVertices);
            for (long j = 0; j < nVertices; j++)
            {
                heapStruct[i]->push(std::numeric_limits<double>::infinity(), j);

                result_SG[(i - 2) * nVertices + j] = -1;
                entering_edge_SG[(i - 2) * nVertices + j] = -1;
            }
        }
        else
        {
            if (i > 0)
            {
                heapStruct[i] = new MinHeap(numNodesOnLayer);
            }
            for (long j = 0; j < numNodesOnLayer; j++)
            {
                if (i == 0)
                {
                    result_CG[j] = -1;
                    predecessors_CG[j] = -1;
                }
                else if (i == 1)
                {
                    heapStruct[i]->push(std::numeric_limits<double>::infinity(), j);
                    result_CG[numNodesOnLayer + j] = -1;
                    predecessors_CG[numNodesOnLayer + j] = -1;
                }
                else
                {
                    heapStruct[i]->push(std::numeric_limits<double>::infinity(), j);
                    result_CG[2 * numNodesOnLayer + j] = -1;
                    predecessors_CG[2 * numNodesOnLayer + j] = -1;
                }
            }
        }
    }

    // run dijkstra
    long nextIndex;
    long level = 0;
    MinHeap *heap = heapStruct[level];
    if (PRUNE && DEBUG_MANIFOLD_DIJKSTRA)
    {
        std::cout << "LineGraph is pruned (no steeringAngles on 3D mesh with more than 90deg allowed)" << std::endl;
    }

    const float progressStepWidth = 1 / 20.0;
    float progressStep = 0;

    

    bool anyEnergyLowerThanUpperBound = false;
    while (!heap->isEmpty())
    {

        if (((float)level / nContour) > progressStep)
        {
            mexPrintf("=");
            // workaround so progress bar actually is visible
            // https://stackoverflow.com/questions/26271154/how-can-i-make-a-mex-function-printf-while-its-running
            mexEvalString("pause(.00001);");
            progressStep += progressStepWidth;
        }

        pair<double, long> current = heap->pop();
        const long currentIdxInLayer = current.second; // This can mean different things based on level! If level =0,1,nContour this is the index of a CGNode, otherwise the index of a product node
        double currentEnergy = current.first;

        
        if (level > 1 && level < nContour)
        {
            result_SG[(level - 2) * nVertices + currentIdxInLayer] = currentEnergy;
        }
        else
        {
            result_CG[getIndexForArrays(level * numNodesOnLayer + currentIdxInLayer, numNodesOnLayer)] = currentEnergy;
        }
        if (currentEnergy <= upperbound)
        {
            if (level > 1 && level < nContour - 1)
            {
                const long current2d = level;
                const long current3d = currentIdxInLayer;

                const long first_edge = entering_edge_SG[current3d + (current2d - 2) * nVertices];

                LGProdNode first_edge_PN = nodesPerLayer[first_edge];

                long previous2d = current2d - 1;
                if (first_edge >= 2 * nEdges3d && first_edge < 4 * nEdges3d)
                {
                    previous2d++;
                }

                // move to next layer on same vertex:
                {
                    long second_edge; 
                    if(first_edge < 2*nEdges3d)
                    {
                        second_edge = first_edge + 4*nEdges3d;
                    }
                    else if(first_edge >= 4*nEdges3d)
                    {
                        second_edge = first_edge;
                    }
                    else
                    {
                        second_edge = first_edge + 2*nEdges3d;
                    }
                    LGProdNode second_edge_PN = nodesPerLayer[second_edge];

                    const long next2d = current2d + 1;
                    updateArrays_SG_SG(first_edge_PN.idxSrc3d, first_edge_PN.idxTrgt3d, second_edge_PN.idxSrc3d, second_edge_PN.idxTrgt3d, previous2d, current2d, next2d,
                                       first_edge, second_edge, currentEnergy, upperbound, entering_edge_SG, numNodesOnLayer,
                                       heapStruct[next2d], result_SG, all_energies, vertices, nVertices, dimVertices, triangles,
                                       nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                }
                for (int k = 0; k <= 1; ++k)
                {
                    for(std::list<std::tuple<long, long>>::iterator neighboring_edge = mesh2.at(current3d)->begin(); neighboring_edge != mesh2.at(current3d)->end(); neighboring_edge++)
                    {
                        const long next2d = current2d + k;
                        long second_edge = std::get<1>(*neighboring_edge);
                        if(k==0)
                        {
                            second_edge += 2*nEdges3d;
                        }
                        LGProdNode second_edge_PN = nodesPerLayer[second_edge];

                        updateArrays_SG_SG(first_edge_PN.idxSrc3d, first_edge_PN.idxTrgt3d, second_edge_PN.idxSrc3d, second_edge_PN.idxTrgt3d, previous2d, current2d, next2d,
                                       first_edge, second_edge, currentEnergy, upperbound, entering_edge_SG, numNodesOnLayer,
                                       heapStruct[next2d], result_SG, all_energies, vertices, nVertices, dimVertices, triangles,
                                       nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                    }
                }
            }
            else if (level < 2)
            {
                anyEnergyLowerThanUpperBound = true;
                LGProdNode currentPN = nodesPerLayer[currentIdxInLayer]; 
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
                const bool isNonDeg = currentIdxInLayer < 2 * nEdges3d;               // diagonal
                const bool isDeg2d = (currentIdxInLayer < 4 * nEdges3d) && !isNonDeg; // horizontal
                const bool isDeg3d = !isDeg2d && !isNonDeg;                           // vertical

                const long next2dSrc = currentPN.idxTrgt2d;
                const long next3dSrc = currentPN.idxTrgt3d;

                const bool isNotLastLayer = next2dSrc <= 1; // "Last" layer means it goes into the level stored as SG, not CG
                const long nextLgStartLayerIdx = getLGStartLayerIndex(next2dSrc, nVertices, nEdges3d);

                // optimized:
                for (auto it : *mesh2.at(next3dSrc))
                {
                    const long next3dTrgt = std::get<0>(it);
                    if (next3dTrgt == currentPN.idxSrc3d)
                    {
                        // no going back i.e. we do not allow going current 3d target be vi if previous source was vi
                        continue;
                    }
                    const long edgeIdx = std::get<1>(it);
                    const long nextLgStartLayerIndexForLoop = getLGStartLayerIndex(next2dSrc, nVertices, nEdges3d);

                    // move to non-deg this or next layer
                    if (isNotLastLayer && (isNonDeg || isDeg3d || isDeg2d))
                    {
                        const long next2dTrgt = next2dSrc + 1;

                        const long nextIdxInLayer = edgeIdx;
                        const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;

                        updateArrays_CG_CG(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                           next2dSrc, next2dTrgt, next3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors_CG, numNodesOnLayer,
                                           heapStruct[next2dSrc], result_CG, all_energies, vertices, nVertices, dimVertices, triangles,
                                           nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                    }

                    // move to deg2d this layer
                    if (isNotLastLayer && isDeg2d)
                    {
                        const long next2dTrgt = next2dSrc;

                        const long nextIdxInLayer = getLGIndexInLayerForDeg2d(edgeIdx, nEdges3d);
                        const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                        updateArrays_CG_CG(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                           next2dSrc, next2dTrgt, next3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors_CG, numNodesOnLayer,
                                           heapStruct[next2dSrc], result_CG, all_energies, vertices, nVertices, dimVertices, triangles,
                                           nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                    }

                    // move to deg2d next layer
                    if (isNotLastLayer && isNonDeg)
                    {
                        const long next2dTrgt = next2dSrc;

                        const long nextIdxInLayer = getLGIndexInLayerForDeg2d(edgeIdx, nEdges3d);
                        const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                        updateArrays_CG_CG(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                           next2dSrc, next2dTrgt, next3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors_CG, numNodesOnLayer,
                                           heapStruct[next2dSrc], result_CG, all_energies, vertices, nVertices, dimVertices, triangles,
                                           nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                    }
                }

                // move to deg3d this or next layer
                if (isNotLastLayer && (isNonDeg || isDeg3d))
                {
                    const long next2dTrgt = next2dSrc + 1;
                    const long next3dTrgt = next3dSrc;
                    const long n3dSrc = currentPN.idxSrc3d;

                    const long nextIdxInLayer = currentIdxInLayer >= 4 * nEdges3d ? currentIdxInLayer : currentIdxInLayer >= 2 * nEdges3d ? currentIdxInLayer + nEdges3dX2
                                                                                                                                          : currentIdxInLayer + nEdges3dX4;
                    const long nextIndex = nextIdxInLayer + nextLgStartLayerIdx;
                    updateArrays_CG_CG(currentIdxInLayer, level, currentPN.idxSrc3d, currentPN.idxTrgt3d, nextIndex, nextIdxInLayer,
                                       next2dSrc, next2dTrgt, n3dSrc, next3dTrgt, currentEnergy, upperbound, predecessors_CG, numNodesOnLayer,
                                       heapStruct[next2dSrc], result_CG, all_energies, vertices, nVertices, dimVertices, triangles,
                                       nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                }

                if (!isNotLastLayer)
                {
                    updateArrays_CG_SG(next3dSrc, currentIdxInLayer, currentEnergy, upperbound, entering_edge_SG, numNodesOnLayer,
                                       heapStruct[next2dSrc], result_SG, all_energies, vertices, nVertices, dimVertices, triangles,
                                       nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                }
            }
            else if (level == nContour - 1) // level=nContour-1, so here we go from SG node i_0 to CG node i_0->i_1.
            {
                const long current2d = level;
                const long current3d = currentIdxInLayer;

                const long first_edge = entering_edge_SG[current3d + (current2d - 2) * nVertices];
                LGProdNode first_edge_PN = nodesPerLayer[first_edge];

                long previous2d = current2d - 1;
                if (first_edge >= 2 * nEdges3d && first_edge < 4 * nEdges3d)
                {
                    previous2d++;
                }

                // move to next layer on same vertex:
                {
                    long second_edge; 
                    if(first_edge < 2*nEdges3d)
                    {
                        second_edge = first_edge + 4*nEdges3d;
                    }
                    else if(first_edge >= 4*nEdges3d)
                    {
                        second_edge = first_edge;
                    }
                    else
                    {
                        second_edge = first_edge + 2*nEdges3d;
                    }

                    LGProdNode second_edge_PN = nodesPerLayer[second_edge];


                    const long next2d = current2d + 1;
                    updateArrays_SG_CG(first_edge_PN.idxSrc3d, first_edge_PN.idxTrgt3d, second_edge_PN.idxSrc3d, second_edge_PN.idxTrgt3d, previous2d, current2d, next2d,
                                       first_edge, second_edge, currentEnergy, upperbound, predecessors_CG, numNodesOnLayer,
                                       heapStruct[level+1], result_CG, all_energies, vertices, nVertices, dimVertices, triangles,
                                       nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                }
                for (int k = 0; k <= 1; ++k)
                {
                    for(std::list<std::tuple<long, long>>::iterator neighboring_edge = mesh2.at(current3d)->begin(); neighboring_edge != mesh2.at(current3d)->end(); neighboring_edge++)
                    {
                        const long next2d = current2d + k;
                        long second_edge = std::get<1>(*neighboring_edge);
                        if(k==0)
                        {
                            second_edge += 2*nEdges3d;
                        }
                        LGProdNode second_edge_PN = nodesPerLayer[second_edge];
                        updateArrays_SG_CG(first_edge_PN.idxSrc3d, first_edge_PN.idxTrgt3d, second_edge_PN.idxSrc3d, second_edge_PN.idxTrgt3d, previous2d, current2d, next2d,
                                           first_edge, second_edge, currentEnergy, upperbound, predecessors_CG, numNodesOnLayer,
                                           heapStruct[level+1], result_CG, all_energies, vertices, nVertices, dimVertices, triangles,
                                           nTriangles, contour, nContour, dimContour, nEdges3d, edges3d, precRotations, landmarks_2d, landmarks_3d, nLandmarks);
                    }
                }
            }
        }

        if (heap->isEmpty() && level < nContour)
        {
            // early stopping since we cannot make progress anymore
            if (!anyEnergyLowerThanUpperBound)
            {
                if (DEBUG_MANIFOLD_DIJKSTRA)
                    std::cout << "Early stopping bc no paths < upperbound found" << std::endl;

                // complete progress bar
                for (int l = level; l < nContour; l++)
                {
                    if (((float)level / nContour) > progressStep)
                    {
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


    for (int i = 0; i < nContour + 1; i++)
    {

        delete heapStruct[i];
    }


    auto t01 = high_resolution_clock::now();
    auto ms_____int = duration_cast<milliseconds>(t01 - t00);
    std::cout << "]  took: " << ms_____int.count() << "ms " << std::endl;
    delete[] heapStruct;
}

#endif /* manifoldDijkstraLG_HPP */
