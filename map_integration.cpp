/*
 *  map_integration.cpp
 *  consensusmap
 *
 *  Created by yonghui on 5/16/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <cassert>
#include "map_integration.h"

namespace consensus_map {

    map_integration::map_integration(individual_mapping_ppl* _in_p_to_array_maps, int _in_number_of_maps)
    {
        assert(_in_p_to_array_maps != NULL);
        number_of_maps = _in_number_of_maps;
        p_to_array_maps = _in_p_to_array_maps;

        srand(time(NULL)); // initialize the random number generator

        // Initialize the lg_clusters object
        cnstrct_lg_clstrs();
        
        // Consolidate the maps into single maps
        lg_clusters.ConsolidateMaps();
        
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void map_integration::cleanup()
    {
        delete[] p_to_array_maps;
        number_of_maps = 0 ;
        p_to_array_maps = NULL;
        return; 
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    map_integration::~map_integration()
    {
        cleanup();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void map_integration::cnstrct_lg_clstrs(){
        // const individual_mapping_ppl* crt_map_ptr;
        int total_number_of_lgs = 0 ; 
        
        for (int ii = 0; ii < number_of_maps; ii++)
        {
            total_number_of_lgs = total_number_of_lgs + p_to_array_maps[ii].get_number_of_lgs();
        }
        
        /*define a temporary vertex type*/
        struct graph_node_str
        {
            const linkage_group_bin* p_to_lg ;
            int ppl_id ;
            int lg_id ;    
        };
        graph_node_str vertices[total_number_of_lgs];
        
        int counter = 0;
        for (int ii = 0 ; ii < number_of_maps ; ii++)
        {
            for (int jj = 0 ; jj < p_to_array_maps[ii].get_number_of_lgs(); jj++)
            {
                vertices[counter].p_to_lg = &((p_to_array_maps[ii].get_lgs())[jj]);
                vertices[counter].ppl_id = ii;
                vertices[counter].lg_id = jj;
                counter = counter + 1;
            }
        }
        
        /*0. construct a boost graph object*/
        typedef adjacency_list <vecS, vecS, undirectedS> Graph;
        Graph lgs_graph(total_number_of_lgs);
        for (int ii = 0 ; ii < total_number_of_lgs; ii++)
        {
            for (int jj = ii+1 ; jj < total_number_of_lgs; jj++ )
            {
                bool overlapping = lgs_intersect(*(vertices[ii].p_to_lg), *(vertices[jj].p_to_lg));
                if (overlapping)
                {
                    add_edge(ii,jj,lgs_graph);
                }
            }
        }

        /*1. Run the connected components algorithm to partition the lgs into clusters*/
        vector<int> cc_ids(total_number_of_lgs, -1);
        int num = connected_components(lgs_graph, &cc_ids[0]);
        
        /*2. Now construct the lg_clusters object*/
        vector<vector<linkage_group_bin*> > lg_groups(num);
        for (int ii = 0 ; ii < total_number_of_lgs; ii++)
        {
            if (cc_ids[ii] >= num)
            {
                cout << "ERROR!, the component id is invalid" << endl;
                assert(cc_ids[ii] >= num); // to crash the program if fail the assert
            }
            lg_groups[cc_ids[ii]].push_back(new linkage_group_bin(*(vertices[ii].p_to_lg)));
        }
        
        vector<LG_CLUSTER*> lgs(num);
        for (int i = 0; i < num; i++) {
            lgs[i] = new LG_CLUSTER();
            lgs[i]->Initialize(lg_groups[i]);
        }
        
        lg_clusters.Initialize(lgs);
        
        // no need to delete those newed objects, since they are taken care of by the lg_clusters object
     };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
