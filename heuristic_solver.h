/*
 *  heuristic_solver.h
 *  consensusmap
 *
 *  Created by yonghui on 1/19/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef HEURISTIC_SOLVER_HEADER
#define HEURISTIC_SOLVER_HEADER

#include <cstdlib>
#include <cassert>
#include <vector>
#include <string>
#include "shortest_path.h"

using namespace std;

namespace consensus_map{
    class H_Solver{
        public:
            // the default constructor
            H_Solver(){
                g_ = NULL;
                n_vertex_ = 0;
            };
            ~H_Solver();
            
            void Initialize(V_Satellite* _g, int _n_vertex) {
                assert((_g != NULL) and (_n_vertex != 0)); // make sure that _g is not null
                g_ = _g;
                n_vertex_ = _n_vertex;
             };
            
            const V_Satellite* get_graph(){
                return g_;
            };
            
            int get_size() {
                return n_vertex_;
            };
            
            // this function returns the id of the marker to be removed
            // If no marker needs to be removed, the function will return a negative number
            int Remove_one_marker();
                                      
        private:

            // a private structure to be used in an internal funciton only
            struct DJ_node {
                int parent;
                bool visited;
                bool reachable;
                double distance;
            };
            
            // initialize the data structure
            void initialize_g();
            void execute_all_pairs_shortest_path();
            // g is set to be null initially 
            // if g is not null, then the destructor is responsible for the destruction of the nodes and edges
            V_Satellite* g_;
            int n_vertex_;            
    };
}

#endif
