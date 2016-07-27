/*
 *  map_integration.h
 *  consensusmap
 *
 *  Created by yonghui on 5/16/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MAP_INTEGRATION_HEADER
#define MAP_INTEGRATION_HEADER

#include "single_population.h"
#include <vector>
#include <set>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <cmath>
#include <queue>
#include <boost/config.hpp>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include "condense_markers_into_bins.h"
#include "ConsensusLG.h"

using namespace std;
using namespace boost;

namespace consensus_map{
    class map_integration
    {
        public:
            
            // _in_p_to_array_maps is supposed to be initialized by the caller. 
            // The callee will take ownership of the _in_p_to_array_maps object after the function call.
            // The object will be destroyed during the destruction of (*this) object.
            // The initializer will do all the consolidation of the individual maps into one single map
            map_integration(individual_mapping_ppl* _in_p_to_array_maps, int _in_number_of_maps);
            ~map_integration();

            void generate_dot_files(){
                lg_clusters.Gen_dot_files();
            };
            
            void Remove_conflicts() {
                lg_clusters.Remove_Conflicts();
            };
            
            void Gen_consensus_dot() {
                lg_clusters.Gen_dot_conflict_free();
            }
            
            void Linearize_Dag() {
                lg_clusters.Linearize_graph();
            }

            void Gen_linear_graph() {
                lg_clusters.Gen_linear_graph();
            }
            
            void Gen_map_chart() { // generate the map in map_chart format
                lg_clusters.Gen_linear_map();
            }
            
            void dump() {
                lg_clusters.Dump();
            };
        // End of public section    
        
        private:
            /*
                The private functions section
            */
           
            // The following function will cluster the set of lg_bins from individual maps into LG clusters
            // The lg_clusters object is populated in the following function
            // The function will be called at the end of the constructor
            void cnstrct_lg_clstrs();
            
            // cleanup will be called in the destructor
            void cleanup();

            /*
                The private data members section
            */
            int number_of_maps; // number of individual populations
            individual_mapping_ppl* p_to_array_maps; // A pointer to an array of individual maps
            
            LG_CLSTR_ENS lg_clusters;
        // End of private section
    };
}//end of namespace
#endif
