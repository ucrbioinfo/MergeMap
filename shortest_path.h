/*
 *  shortest_path.h
 *  consensusmap
 *
 *  Created by yonghui on 11/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef SHORTEST_PATH_HEADER
#define SHORTEST_PATH_HEADER

#include <cstdlib>
#include <cassert>
#include <vector>
#include <string>
using namespace std;

namespace consensus_map{

    // The satellite data associated with an edge
    struct E_Satellite {
        int to_id;
        bool special;
    };
    
    // The satellite data associate with an node
    struct V_Satellite {
        double weight; // the cost of the node
        double dual_val; // dual variable associated with the node
        double total_flow; // total flow through the node 
        double nor_dual_val; // the dual variable normalized to [0,1]
        bool to_delete;
        string name; // for the purpose of debugging
        int num_edges;
        E_Satellite* edges;
    };
    
    class CG {
        public:
            CG(){
                // initialize g to be null
                g_ = NULL;
                n_vertex_ = 0;
            };
            
            ~CG();
            
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
            
            void solve_lp(double epsilon);
            
            void dump();
            
            V_Satellite* construct_sub_graph(const vector<int>& vertices) const;
            
        private:
        
            // a private structure to be used in an internal funciton only
            struct DJ_node {
                int parent;
                bool visited;
                bool reachable;
                double distance;
            };
            
            // find the shortest path and return it as a vector
            void find_shortest_path(double& crt_cost, vector<int>& opt_path);
            
            // compute dual values
            void compute_dual_vals(double alpha);
            
            // compute lambda
            double compute_lambda();
            
            // update primal solution
            void update_solution(double delta, vector<int> & path);

            // initialize rounding step
            void initialize_lp_round();
            
            // round the fractional solution to integral solutions
            void round_lp_randomized(vector<int> nodes);            
            void round_lp_greedy(vector<int> nodes);
            // will return true if there exists none-special edge in the sub-graph induced by nodes
            bool exist_no_special_edge(const vector<int>& nodes);
            bool valid_solution();
            // minimize the existing solution to a minimal solution by using some hill-climbing approach
            void minimize_sol(); 

            void round_lp();
           
            // g is set to be null initially 
            // if g is not null, then the destructor is responsible for the destruction of the nodes and edges
            V_Satellite* g_;
            int n_vertex_;
    };
}
#endif //SHORTEST_PATH_HEADER
