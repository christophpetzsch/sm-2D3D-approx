/**
		@author Zorah Lähner (laehner@in.tum.de)
    @version 03.2016
*/
#ifndef buildMeshLG_HPP
#define buildMeshLG_HPP

#include <vector>
#include <list>
#include <set>
#include <tuple>
#include <cassert>

using std::vector;
using std::list;


struct EDGE {
    int midx0; 
    int midx1; 
    int mmaxNumEdges;
    EDGE(int idx0, int idx1, int maxNumEdges) {
        midx0 = idx0; 
        midx1 = idx1; 
        mmaxNumEdges = maxNumEdges; 
    }
    
    bool operator< (const EDGE& otherEdge) const {
        return (this->midx0 * this->mmaxNumEdges + this->midx1) < 
                (otherEdge.midx0 * otherEdge.mmaxNumEdges + otherEdge.midx1);
    }
};

/* creates the adjacency list of triangles into mesh
 * vector should have appropriate size and be filled with initialised lists already
 *
 * List has the following format: 
 * vector[vertexId] = list of connected vertexIds on triangle mesh
 */
void buildMesh(const double* triangles, int nTriangles, vector<list<long>* >& mesh, 
        vector<list<std::tuple<long, long>>* >& mesh2, 
        vector<std::tuple<long, long>>& edges3d) {
	// add all edges
	for(int i=0; i < nTriangles; i++) {
		for (int j=0; j < 3; j++) {
			mesh.at(int(triangles[3*i+j]))->push_back(int(triangles[3*i+((j + 1) % 3)]));
			mesh.at(int(triangles[3*i+j]))->push_back(int(triangles[3*i+((j + 2) % 3)]));
		}
	}
	// remove double entries
	for(int i=0; i < mesh.size(); i++) {
		mesh.at(i)->sort();
		mesh.at(i)->unique();
	}  
    
    
    std::set<EDGE> edges;
    //edges.reserve(6 * mesh.size());
    const int maxNumEdges = 3 * nTriangles; 
    
	for(int i=0; i < nTriangles; i++) {
        for (int j = 0; j < 3; j++) {
            const int idx0 = int(triangles[3*i + ((j + 0) % 3)]);
            const int idx1 = int(triangles[3*i + ((j + 1) % 3)]);
            
            EDGE firstDir = EDGE(idx0, idx1, maxNumEdges);
            EDGE secondDir = EDGE(idx1, idx0, maxNumEdges);
           
            if (edges.find(firstDir) != edges.end()) {
                // edge already in set
                continue;
            }
            else if (edges.find(secondDir) != edges.end()) {
                // other edge direction already in set
                continue;
            }
            else {
                edges.insert(firstDir);
            }
        }
	}
    
    int idx = 0;
    const int nEdges = edges.size();
    edges3d.reserve(nEdges);
    for(const EDGE &edge : edges) {        
        mesh2.at(edge.midx0)->push_back(std::tuple<long, long>{edge.midx1, idx});
        mesh2.at(edge.midx1)->push_back(std::tuple<long, long>{edge.midx0, nEdges + idx});
        idx++;
        
        edges3d.push_back(std::tuple<long, long>{edge.midx0, edge.midx1});
    }
    
    /*for(int i=0; i < mesh2.size(); i++) {
        for (auto it : *mesh2.at(i))
            std::cout << std::get<1>(it) << ": " << i << " -> " << std::get<0>(it) << "; ";
                
        std::cout << std::endl;
	} */

}


#endif /* buildMeshLG_HPP */
