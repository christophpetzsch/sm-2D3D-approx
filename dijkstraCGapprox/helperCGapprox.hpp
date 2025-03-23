/**
    Helper functions..
    @author Zorah Lähner (laehner@in.tum.de)
    @version 03.2016
*/
#include <cassert>

struct LGProdNode
{
    long idxSrc2d;
    long idxTrgt2d;
    long idxSrc3d;
    long idxTrgt3d;
    LGProdNode(long iidxSrc2d, long iidxTrgt2d, long iidxSrc3d, long iidxTrgt3d)
    {
        idxSrc2d = iidxSrc2d;
        idxTrgt2d = iidxTrgt2d;
        idxSrc3d = iidxSrc3d;
        idxTrgt3d = iidxTrgt3d;
    }
    LGProdNode()
    {
    }
};


long getVertex(long index, long nVertices)
{
    return long(index % nVertices);
}

long getContour(long index, long nVertices)
{
    return index / nVertices;
}

long getIndex(long vertex, long contour, long nVertices)
{
    return contour * nVertices + vertex;
}

/* use mesh2 vectorListTuple to get edge idx
 * - each vector element belongs to a 3d vertex
 * - each vector element is a list of connected 3d vertices to respective 3d vertex
 * - each element in the list is a tuple containing first the idx of adjacent 3d vertex
 *   and second the idx in the list of edges
 */
long getEdgeIdx(long idxSrc3d, long idxTrgt3d, const std::vector<std::list<std::tuple<long, long>> *> mesh2)
{

    for (auto const &it : *mesh2.at(idxSrc3d))
    {
        if (std::get<0>(it) == idxTrgt3d)
        {
            return std::get<1>(it);
        }
    }
    // we should never reach this code
    assert(false);
    std::cout << "Reached code section which should not be reached" << std::endl;
    return 0;
}

std::tuple<long, long> getEdge(long edgeIdx, const vector<std::tuple<long, long>> edges3d)
{
    const int nEdges = edges3d.size();
    bool swapSrcTrgt = false;
    if (edgeIdx >= nEdges)
    {
        edgeIdx -= nEdges;
        swapSrcTrgt = true;
    }
    const std::tuple<long, long> edge = edges3d.at(edgeIdx);
    const long idxSrc = std::get<0>(edge);
    const long idxTrgt = std::get<1>(edge);

    if (swapSrcTrgt)
    {
        return std::tuple<long, long>{idxTrgt, idxSrc};
    }
    else
    {
        return std::tuple<long, long>{idxSrc, idxTrgt};
    }
}

/* Order of each line graph layer:
 * 1) non-degenerate edges in first 3d edge orientation
 * 2) non-degenerate edges in second 3d edge orientation
 * 3) degenerate 2d edges in first 3d edge orientation
 * 4) degenerate 2d edges in second 3d edge orientation
 * 5) degenerate 3d edges
 */
inline long getLGIndexInLayerForDeg2d(const long edgeIdx, const long nEdges3d)
{
    return edgeIdx + 2 * nEdges3d;
}

// inline long getLGIndexInLayerForDeg3d(const long idxSrc3d, const long nEdges3d) {
//     return idxSrc3d + 4 * nEdges3d;
inline long getLGIndexInLayerForDeg3d(const long edgeIdx, const long nEdges3d)
{
    return edgeIdx + 4 * nEdges3d;
}

inline long getLGIndexInLayerForNonDeg(const long edgeIdx)
{
    return edgeIdx;
}

/*inline long getLGIndexInLayerProvidingEdgeIdx(long edgeIdx, long idxSrc2d, long idxTrgt2d, long idxSrc3d, long idxTrgt3d, const long nVertices3d, const long nEdges3d,
        const std::vector<std::list<std::tuple<long, long>>*> mesh2) {
    int idx = 0;
    // 3 & 4)
    if (idxSrc2d == idxTrgt2d) {
        idx = getLGIndexInLayerForDeg2d(edgeIdx, nEdges3d);
    }
    // 5
    else if (idxSrc3d == idxTrgt3d) { // ???
        //idx = getLGIndexInLayerForDeg3d(idxSrc3d, nEdges3d);
        idx = getLGIndexInLayerForDeg3d(edgeIdx, nEdges3d);
    }
    // 1 & 2)
    else {
        idx = getLGIndexInLayerForNonDeg(edgeIdx);
    }
    return idx;
}*/

/*long getLGIndexInLayer(long idxSrc2d, long idxTrgt2d, long idxSrc3d, long idxTrgt3d, const long nVertices3d, const long nEdges3d,
        const std::vector<std::list<std::tuple<long, long>>*> mesh2) {
    const int edgeIdx = idxSrc3d == idxTrgt3d ? -1 : getEdgeIdx(idxSrc3d, idxTrgt3d, mesh2);
    return getLGIndexInLayerProvidingEdgeIdx(edgeIdx, idxSrc2d,idxTrgt2d, idxSrc3d, idxTrgt3d, nVertices3d, nEdges3d, mesh2);
}*/

long getLGStartLayerIndex(long idx2d, const long nVertices3d, const long nEdges3d)
{
    // return idx2d * (2 * nEdges3d + nVertices3d + 2 * nEdges3d);
    return idx2d * (6 * nEdges3d);
}

/*long getLGIndex(long idxSrc2d, long idxTrgt2d, long idxSrc3d, long idxTrgt3d, const long nVertices3d, const long nEdges3d,
        const std::vector<std::list<std::tuple<long, long>>*> mesh2) {

    int idx = getLGStartLayerIndex(idxSrc2d, nVertices3d, nEdges3d);
    // needs change
    return idx + getLGIndexInLayer(idxSrc2d, idxTrgt2d, idxSrc3d, idxTrgt3d, nVertices3d, nEdges3d, mesh2);
}*/

LGProdNode getProductNodeInLayer(long lgIdxInLayer, long idxSrc2d,
                                 const vector<std::tuple<long, long>> edges3d, const int nVertices3d)
{

    const int nEdges3d = edges3d.size();
    LGProdNode lgpn;
    lgpn.idxSrc2d = idxSrc2d;
    // 1 & 2)
    if (lgIdxInLayer < 2 * nEdges3d)
    {
        lgpn.idxTrgt2d = lgpn.idxSrc2d + 1;
        std::tuple<long, long> edge = getEdge(lgIdxInLayer, edges3d);
        lgpn.idxSrc3d = std::get<0>(edge);
        lgpn.idxTrgt3d = std::get<1>(edge);
    }
    // 3 & 4)
    else if (lgIdxInLayer < 4 * nEdges3d)
    {
        lgIdxInLayer -= 2 * nEdges3d;
        lgpn.idxTrgt2d = lgpn.idxSrc2d;
        std::tuple<long, long> edge = getEdge(lgIdxInLayer, edges3d);
        lgpn.idxSrc3d = std::get<0>(edge);
        lgpn.idxTrgt3d = std::get<1>(edge);
    }
    // 5)
    else
    {
        lgIdxInLayer -= 4 * nEdges3d;
        lgpn.idxTrgt2d = lgpn.idxSrc2d + 1;
        // lgpn.idxSrc3d = lgIdxInLayer;
        // lgpn.idxTrgt3d = lgpn.idxSrc3d;
        std::tuple<long, long> edge = getEdge(lgIdxInLayer, edges3d);
        lgpn.idxSrc3d = std::get<0>(edge);
        lgpn.idxTrgt3d = std::get<1>(edge);
    }

    return lgpn;
}

LGProdNode getProductNode(long lgIdx,
                          const vector<std::tuple<long, long>> edges3d, const int nVertices3d)
{

    const int nEdges3d = edges3d.size();
    // const int layer = lgIdx / (2 * nEdges3d + nVertices3d + 2 * nEdges3d);
    const int layer = lgIdx / (6 * nEdges3d);
    // return getProductNodeInLayer(lgIdx - layer * (2 * nEdges3d + nVertices3d + 2 * nEdges3d), layer, edges3d, nVertices3d);
    return getProductNodeInLayer(lgIdx - layer * (6 * nEdges3d), layer, edges3d, nVertices3d);
}

long getLGSrc3dVertex(int lgNodeIdx,
                      const vector<std::tuple<long, long>> edges3d, const int nVertices3d)
{
    LGProdNode lgpn = getProductNode(lgNodeIdx, edges3d, nVertices3d);
    if (lgNodeIdx - getLGStartLayerIndex(lgpn.idxSrc2d, nVertices3d, edges3d.size()) >= 4 * edges3d.size())
        return lgpn.idxTrgt3d;
    return lgpn.idxSrc3d;
}

long getLGTrgt3dVertex(int lgNodeIdx,
                       const vector<std::tuple<long, long>> edges3d, const int nVertices3d)
{
    LGProdNode lgpn = getProductNode(lgNodeIdx, edges3d, nVertices3d);
    return lgpn.idxTrgt3d;
}
/*void testLGHelperFunctions(const std::vector<std::list<std::tuple<long, long>>*> mesh2,
                           const vector<std::tuple<long,long>> edges3d,
                           int nVertices3d,
                           int nVertices2d) {

    const int nEdges3d = edges3d.size();
    std::cout << "nEdges3d = " << nEdges3d << std::endl;
    std::cout << "Testing edgeidx to edge conversion..." << std::endl;
    for (int e = 0; e < 2 * nEdges3d; e++) {
        std::tuple<long, long> edge = getEdge(e, edges3d);
        const int idxSrc3d  = std::get<0>(edge);
        const int idxTrgt3d = std::get<1>(edge);
        int e_ = getEdgeIdx(idxSrc3d, idxTrgt3d, mesh2);
        if (e != e_) {
            std::cout << "ERROR in edgeIdx to edge conversion: idx = " << e << " e_ = " << e_ << std::endl;
            return;
        }
    }
    std::cout << "Done." << std::endl;

    std::cout << "Testing LG product nodes to product node conversion..." << std::endl;
    for (int layer = 0; layer < nVertices2d; layer++) {
        for (int i = 0; i < 4 * nEdges3d + nVertices3d; i++) {
            int idx = layer * (4 * nEdges3d + nVertices3d) + i;
            LGProdNode lgpn = getProductNode(idx, edges3d, nVertices3d);
            int i_ = getLGIndex(lgpn.idxSrc2d, lgpn.idxTrgt2d, lgpn.idxSrc3d, lgpn.idxTrgt3d, nVertices3d, nEdges3d, mesh2);
            if (idx != i_) {
                std::cout << "ERROR in productNodeIdx to productNode conversion: idx = " << idx << " i_ = " << i_ << std::endl;
                return;
            }
        }
    }
    std::cout << "Done." << std::endl;
}*/

double curvature(double p1, double p2, double p3, double q1, double q2, double q3, double r1, double r2, double r3)
{
    double a = sqrt((p1 - q1) * (p1 - q1) + (p2 - q2) * (p2 - q2) + (p3 - q3) * (p3 - q3));
    double b = sqrt((q1 - r1) * (q1 - r1) + (q2 - r2) * (q2 - r2) + (q3 - r3) * (q3 - r3));
    double c = sqrt((r1 - p1) * (r1 - p1) + (r2 - p2) * (r2 - p2) + (r3 - p3) * (r3 - p3));

    double p = (a + b + c) / 2;
    double k = sqrt(p * (p - a) * (p - b) * (p - c));

    return (4 * k) / (a * b * c);
}

long getLocalIndex(long vertex, long nVertices)
{
    return vertex % nVertices;
}

long getIndexForArrays(long global_edge, const long numNodesOnLayer)
{
    if (global_edge < 2 * numNodesOnLayer)
    {
        return global_edge;
    }
    else
    {
        return (global_edge % numNodesOnLayer + 2 * numNodesOnLayer);
    }
}

void calcPath(long end, const long *predecessors_CG, const long *entering_edge_SG, double *result_CG, double *result_SG, list<long> &path, const long nContour, const int nVertices3d, const int nEdges3d, const vector<std::tuple<long, long>> &edges3d, const LGProdNode *nodesPerLayer)
{
    long local_edge = end-12*nEdges3d;

    path.push_front(local_edge + nContour*6*nEdges3d);
    local_edge = predecessors_CG[local_edge+12*nEdges3d];

    LGProdNode local_edge_PN = nodesPerLayer[local_edge];


    for (long level = nContour - 2; level > 1; --level)
    {
        path.push_front(local_edge + level*6*nEdges3d);

        if (local_edge >= 2 * nEdges3d && local_edge < 4 * nEdges3d)
        {
            long local_vertex = local_edge_PN.idxSrc3d;
            local_edge = entering_edge_SG[(level-2) * nVertices3d + local_vertex];
        }
        else if(local_edge < 2*nEdges3d)
        {
            long local_vertex = local_edge_PN.idxSrc3d;
            local_edge = entering_edge_SG[(level-2) * nVertices3d + local_vertex];
        }
        else 
        {
            long local_vertex = local_edge_PN.idxTrgt3d;
            local_edge = entering_edge_SG[(level-2) * nVertices3d + local_vertex];
        }
        if(local_edge >= 2 * nEdges3d && local_edge < 4 * nEdges3d) ++level; // we will decrease level soon even though we took a horizontal edge, thus we increase it now by one
        local_edge_PN = nodesPerLayer[local_edge];
    }

    path.push_front(local_edge + 6*nEdges3d);
    local_edge = predecessors_CG[local_edge+6*nEdges3d];

    while(local_edge >= 2 * nEdges3d && local_edge < 4 * nEdges3d)
    {
        path.push_front(local_edge + 6*nEdges3d);
        local_edge = predecessors_CG[local_edge+6*nEdges3d];
    }

    path.push_front(local_edge);
    local_edge = predecessors_CG[local_edge];

    while (local_edge > -1)
    {
        path.push_front(local_edge);
        local_edge = predecessors_CG[local_edge];
    }
}
