/*
 *  heuristic_solver.cpp
 *  consensusmap
 *
 *  Created by yonghui on 1/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include "heuristic_solver.h"
#include <iostream>

namespace consensus_map {
    H_Solver::~H_Solver(){
        if (g_ != NULL) {
            assert(n_vertex_ != 0);
            for (int ii = 0; ii < n_vertex_; ii++) {
                delete[] g_[ii].edges;
            }
            delete[] g_;
            g_ = NULL;
            n_vertex_ = 0;
        }
    } // end of CG::~CG()
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void H_Solver::initialize_g(){
        for (int ii = 0; ii < n_vertex_; ii++) {
            g_[ii].to_delete = false;
            g_[ii].dual_val = 1.0; // the shortest paths are un-weighted
            g_[ii].total_flow = 0.0;
            g_[ii].nor_dual_val = 0.0;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int H_Solver::Remove_one_marker(){
        // initialize the graph
        initialize_g();
        // find all the shortest cycles 
        execute_all_pairs_shortest_path();
        
        vector<int> possible_markers;
        for (int ii = 0; ii < n_vertex_; ii++) {
            if (g_[ii].total_flow > 0) {
                possible_markers.push_back(ii);
            }
        }
        
        // it is expected that at least one marker should be removed, otherwise, the routine should not be called at all
        if(possible_markers.size() == 0){
            cout << "nothing to be deleted" << endl;
            return -1; 
        }; 
        int to_remove = -1;
        double max_satuation = numeric_limits<double>::max();
        for (int jj = 0; jj < possible_markers.size(); jj++) {
            int node_id = possible_markers[jj];
            if ((g_[node_id].weight / g_[node_id].total_flow) < max_satuation ) {
                to_remove = node_id;
                max_satuation = g_[node_id].weight / g_[node_id].total_flow; 
            }
        }
        assert(to_remove != -1);
        return to_remove;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    void H_Solver::execute_all_pairs_shortest_path(){
        
        DJ_node* node_status = new DJ_node[n_vertex_];
        // run dijastra's algorithm starting from every node
        for (int ii = 0; ii < n_vertex_; ii++) {
            // initialize the node_status array
            for (int jj = 0; jj < n_vertex_; jj++) {
                node_status[jj].parent = -1;
                node_status[jj].visited = false;
                node_status[jj].reachable = false;
                node_status[jj].distance = numeric_limits<double>::max();
            }
            // special vertices
            // it is expected that the sinks set is non-empty
            vector<int> sinks;
            // initialize the first node;
            node_status[ii].visited = true;
            node_status[ii].reachable = true;
            node_status[ii].distance = g_[ii].dual_val;
            for (int jj = 0; jj < g_[ii].num_edges; jj++) {
                int dest_id = g_[ii].edges[jj].to_id;
                if (g_[ii].edges[jj].special) {
                    sinks.push_back(dest_id);
                } else {
                    node_status[dest_id].distance = node_status[ii].distance + g_[dest_id].dual_val;
                    node_status[dest_id].parent = ii;
                    node_status[dest_id].reachable = true;
                }
            }
            
            
            while (true) {
                // find the closest node that still has not been visited
                double min_dist = numeric_limits<double>::max();
                int dest_id = -1;
                for (int kk = 0; kk < n_vertex_; kk++) {
                    if ((not node_status[kk].visited) and 
                        (node_status[kk].reachable) and 
                        (node_status[kk].distance < min_dist)) {
                        min_dist = node_status[kk].distance;
                        dest_id = kk;
                    }
                }
                if (dest_id != -1) {
                    node_status[dest_id].visited = true;
                    for (int ll = 0; ll < g_[dest_id].num_edges; ll++) {
                        int ll_id = g_[dest_id].edges[ll].to_id;
                        if (not node_status[ll_id].visited) {
                            if (node_status[ll_id].distance > (node_status[dest_id].distance + g_[ll_id].dual_val)) {
                                node_status[ll_id].distance = node_status[dest_id].distance + g_[ll_id].dual_val;
                                node_status[ll_id].parent = dest_id;
                                node_status[ll_id].reachable = true;
                            }
                        }
                    }
                } else {
                    break;
                }
            } // end while
            
            for (vector<int>::iterator sink = sinks.begin();
                 sink != sinks.end();
                 ++sink) {
                if (node_status[*sink].visited) {
                    g_[*sink].total_flow = g_[*sink].total_flow + 1; 
                    int parent = node_status[*sink].parent;
                    while (parent > -1) {
                        g_[parent].total_flow = g_[parent].total_flow + 1; 
                        parent = node_status[parent].parent;
                    }                     
                }
            }
                       
        } // end of for (dijastra to find shortest path from each node)
        
        delete[] node_status;
        return;        
    };
}
