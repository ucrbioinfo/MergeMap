/*
 *  ConsensusLG.h
 *  consensusmap
 *
 *  Created by yonghui on 10/15/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONSENSUSLG_HEADER
#define CONSENSUSLG_HEADER

#include "single_population.h"
#include "condense_markers_into_bins.h"
#include "fragment.h"
#include "SCC_Drawer.h"

namespace consensus_map{
    
    class SCC_Drawer;
    class linkage_group_graph_representation;
    class linkage_group_bin;
    // The following class is intended to encapsulate related LGs from individual maps that 
    // should be consolidated into one single consensu Linkage Group.
    // raw_lgs_ are supposed to be in the right order already.
    // Each _raw_lgs[i] comes from one population.
    class LG_CLUSTER{
        public:
            // Define the constructor and the destructor
            LG_CLUSTER();
            ~LG_CLUSTER();
            // The following initializer will copy the elements from the _raw_lgs vector to 
            // the raw_lgs_ vector and will then set those pointers in _raw_lgs to NULL.
            // In other words, the ownership of those linkage_group_bin objects will be transfered to the (this) object
            // to be constructed.
            void Initialize(vector<linkage_group_bin*>& _raw_lgs);

            
            // Several simple accessors
            int num_raw_lgs(){
                return raw_lgs_.size();
            };
            
            // The following function will remove conflicting nodes from individual maps 
            // and will re-build the maps
            void remove_conflicts();
                                    
            // Consolidate the LG cluster into one single LG in graph representation
            void ConsolidateCluster(){
                SimplifyCluster();
                gen_consensus_LG();
            };
            
            // Dump some statistics
            void Dump();
            
            // Generate a fragment from each LG according to a strongly connected component
            FragClstr* GenFragClstr(const vector<string>& scc);
            
            // Generate a fragment drawer
            SCC_Drawer* GenDrawer(int scc_id, bool direction, vector<string>& scc);
            
            // Generate the Dot files 
            void Gen_dot_files(string file_path);
            
            // accumulate conflicts from different fragments
            void accumulate_conflicts(const vector<pair<string, int> >& conflicts);
            
            void Linearize_Dag();
            
            void Gen_linear_map(string file_path);
            
            string Gen_map_chart(); // generate the linkage group in mapchart format
            
        private:
            // Simplify the LG_CLUSTER. 
            // In this step, markers are condensed into supernodes to simplify the problem of building
            // consensus map and will make the final DAG much cleaner.
            void SimplifyCluster(){
                fix_orientation();
                condense_into_sn_scc();
                condense_into_sn();
            };
            
            // The following three functions are to be called sequentially in order to simplify an LG_CLUSTER
            // Step 1: figure out the right orientation of each individual linkage_group_bin,
            //         and flip some of them if necessary. fix_orientation_help is a helper for the fix_orientation 
            //         function and is to be called by fix_orientation
            void fix_orientation();
            vector<bool> fix_orientation_help(const vector<vector<double> >& pair_wise_tau);
            // Step 2: Condense the markers into the super node representation. 
            //         Two nodes are condensed into one super-node if they cosegerate in some ppl,
            //         and there is no information to break them. 
            void condense_into_sn();
            // Step 3: Condense the markers from the same Strongly Connected Component into supernode.
            //         This step is useful when a SCC involves a lot of markers. Clearly the previous step will
            //         never condense any markers since there are conflicting relationships between them. 
            //         However, it is possible that there are pairs of markers which always co-appear in 
            //         the same bin
            void condense_into_sn_scc();
            
            // This step will actually generate the consensus LG in adjacency list representation
            // The resulting graph will be pointed to by consensus_LG
            // This function can be called multiple times
            void gen_consensus_LG();
            
            
            // some convenient functions 
            bool IdenticalSig(const vector<int>& sig1, const vector<int>& sig2){
                assert(sig1.size() == sig2.size());
                for (int ii = 0; ii < sig1.size(); ii++) {
                    if (sig1[ii] != sig2[ii]) {
                        return false;
                    }
                }
                return true;
            }
            
            
            // those linkage_group_bin objects pointed by the elements in the vector are actually owned by 
            // this class. They should be freed when the object is destructed, otherwise, 
            // memory leak will result. 
            // This vector is initialized to NULLs. 
            vector<linkage_group_bin*> raw_lgs_;
            
            // data generated after consolidation of the LG cluster
            linkage_group_graph_representation* consensus_LG;
            
            // an internal state variable
            bool conflicts_computed;
            vector<vector<string> > nodes_to_remove; 
    };
    
    // The following class represents all the LG_CLUSTER to be assembled 
    class LG_CLSTR_ENS{
        public:
            // Define the constructor and the destructor
            LG_CLSTR_ENS();
            ~LG_CLSTR_ENS();
            // The following function will copy the elements from the _lgs vector to the lgs_ vector.
            // The elements in the _lgs vector will be set to NULL when the function finishes.
            // In other words, the initializer will transfer the ownship of LG's from _lgs to (this) object
            void Initialize(vector<LG_CLUSTER*>& _lgs);
            
            // Define several accessor
            int num_clusters() const{
                return lg_clusters_.size();
            };
            
            void ConsolidateMaps(){
                for (int i = 0; i < lg_clusters_.size(); i++) {
                    (lg_clusters_[i])->ConsolidateCluster();
                }                
            };
            
            void Remove_Conflicts(){
                for (int i = 0; i < lg_clusters_.size(); i++) {
                    (lg_clusters_[i])->remove_conflicts();
                }                
            }
            
            void Dump(){
                for (int i = 0; i < lg_clusters_.size(); i++) {
                    (lg_clusters_[i])->Dump();
                }
            };
            
            void Gen_dot_files();
            
            void Linearize_graph(){
                for (int ii = 0; ii < lg_clusters_.size(); ii++) {
                    lg_clusters_[ii]->Linearize_Dag();
                }
            }
            
            void Gen_dot_conflict_free();
            
            void Gen_linear_graph(); // output the lgs in dot file format
            
            void Gen_linear_map(); // output the lgs in the mapchart format
                        
        private:
            vector<LG_CLUSTER*> lg_clusters_;
    };
}

#endif
