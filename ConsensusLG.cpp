/*
 *  ConsensusLG.cpp
 *  consensusmap
 *
 *  Created by yonghui on 10/15/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <cassert>
#include "ConsensusLG.h"

namespace consensus_map {

    LG_CLUSTER::LG_CLUSTER(){
        consensus_LG = NULL;
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLUSTER::Initialize(vector<linkage_group_bin*>& _raw_lgs) {
        for (int i = 0; i < _raw_lgs.size(); i++) {
            // it is assumed that none of the pointers is NULL
            assert(_raw_lgs[i] != NULL);
            raw_lgs_.push_back(_raw_lgs[i]);
            _raw_lgs[i] = NULL;
        }
        conflicts_computed = false;
        nodes_to_remove.resize(_raw_lgs.size());
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    LG_CLUSTER::~LG_CLUSTER() {
        for (int i = 0; i < raw_lgs_.size(); i++) {
            assert(raw_lgs_[i] != NULL);
            delete raw_lgs_[i];
        };
        if (consensus_LG != NULL) {
            delete consensus_LG;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLUSTER::fix_orientation() {
        int lg_size = raw_lgs_.size();
        vector<vector<double> > pair_wise_tau_matrix(lg_size, vector<double>(lg_size));
        vector<vector<double> > pair_wise_tau_matrix_rev(lg_size, vector<double>(lg_size));
        vector<vector<int> > pair_wise_common_markers(lg_size, vector<int>(lg_size));
        for (int mm = 0 ; mm < lg_size; mm++) {
            for (int nn = 0 ; nn < lg_size; nn++) {
                pair<double,double> tau = pairwise_tau(*(raw_lgs_[mm]),
                                                       *(raw_lgs_[nn]),
                                                       pair_wise_common_markers[mm][nn]);
                pair_wise_tau_matrix[mm][nn] = tau.first;
                pair_wise_tau_matrix_rev[mm][nn] = tau.second;
            }
        }

        vector<bool> orientations = fix_orientation_help(pair_wise_tau_matrix);
        for (int ll= 0 ; ll < lg_size; ll++) {
            if (orientations[ll] == false) {
                // flip the orientation of an LG if necessary
                raw_lgs_[ll]->flip_orientation();
            }
        }

        // Print out some debugging information
        cout << "pair-wise tau for consensus lg:" << endl;
        char buffer[100];
        for (int aa = 0; aa < lg_size; aa++ ) {
            for (int bb = 0 ; bb < lg_size; bb++) {
                cout << '\t';
                sprintf(buffer, "%+.3f", pair_wise_tau_matrix[aa][bb]);
                cout << buffer;
                cout << ",";
                sprintf(buffer, "%+.3f", pair_wise_tau_matrix_rev[aa][bb]);
                cout << buffer;
                sprintf(buffer, "%3d", pair_wise_common_markers[aa][bb]);
                cout << '(' << buffer << ')';
            }
            cout << endl;
        }
        cout << endl;
        // End printing debugging information
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    vector<bool> LG_CLUSTER::fix_orientation_help(const vector<vector<double> >& pair_wise_tau)
    {
        int number_of_lgs = pair_wise_tau.size();
        vector<bool> visitted(number_of_lgs, false);
        vector<int> orientations(number_of_lgs, 1);
        vector<bool> to_return(number_of_lgs, true);
        std::queue<int> nodes_to_be_visitted;
        for (int ii = 0 ; ii < number_of_lgs; ii++)
        {
            if (visitted[ii] == false)
            {
                visitted[ii] = true;
                nodes_to_be_visitted.push(ii); // as long as the node gets in the queue, it is visitted
            }
            while (not nodes_to_be_visitted.empty())
            {
                int next_id = nodes_to_be_visitted.front();
                nodes_to_be_visitted.pop();
                for (int jj = 0 ; jj < number_of_lgs; jj++)
                {
                    if((not visitted[jj]) and (abs(pair_wise_tau[next_id][jj]) > EPSILON))
                    {
                        visitted[jj] = true;
                        if (pair_wise_tau[next_id][jj] < 0)
                        {
                            orientations[jj] = -orientations[next_id];
                        }
                        else 
                        {
                            orientations[jj] = orientations[next_id];
                        }
                        nodes_to_be_visitted.push(jj);
                    }
                }
            }
        }
        /*check to make sure that the orientation is consistent*/
        for (int ii = 0 ; ii < number_of_lgs; ii++)
        {
            for (int jj = 0 ; jj < number_of_lgs; jj++)
            {
                if (pair_wise_tau[ii][jj] * orientations[ii] * orientations[jj] < -2*EPSILON)
                {
                    cout << "ERROR! no consistent orientation can be found" << endl;
                }
            }
        }
        
        for (int ii = 0; ii < number_of_lgs; ii++)
        {
            if (orientations[ii] < 0)
            {
                to_return[ii] = false;
            }
        }
        return to_return;
        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLUSTER::condense_into_sn() {
        /*
            It is assumed that the orientation for the individual LGs has already been fixed
        */
        int num_lgs = raw_lgs_.size();
        
        // The consensus_lgs_dag needs to be deleted at the end
        linkage_group_DAG* consensus_lgs_dag = NULL;
        
        // Step 0: construct the linkage_group_DAG object
        assert(raw_lgs_.size() > 0);
        consensus_lgs_dag = raw_lgs_[0]->construct_linkage_group_DAG(true);

        for (int mm = 1; mm < num_lgs; mm++) {
            linkage_group_DAG* tmp_lg_dag1 = raw_lgs_[mm]->construct_linkage_group_DAG(true);
            linkage_group_DAG* tmp_lg_dag2 = consensus_lgs_dag;            
            consensus_lgs_dag = combine_lg_DAGs(tmp_lg_dag1, tmp_lg_dag2);
            delete tmp_lg_dag1;
            delete tmp_lg_dag2;
        }

        vector<const linkage_group_bin*> in_bins(num_lgs);
        for (int ii2 = 0; ii2 < num_lgs; ii2++)
        {
            in_bins[ii2] = raw_lgs_[ii2];
        }

        // Step 1: consolidate the markers into super nodes
        vector<linkage_group_bin*> out_bins(num_lgs);
        consolidate_markers_into_bins(in_bins, num_lgs, out_bins, *consensus_lgs_dag);
        for (int ii2 = 0; ii2 < num_lgs; ii2++)
        {
            delete raw_lgs_[ii2];
            raw_lgs_[ii2] = out_bins[ii2];
        }
        
        // Deallocation of the memory at the end of the function
        assert(consensus_lgs_dag != NULL);
        delete consensus_lgs_dag;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLUSTER::condense_into_sn_scc() {

        // To be conservative, two markers are gonna be condensed into one supernode if they 
        // always co-appear in the same bin in all the populations
        int num_lgs = raw_lgs_.size();
        
        // step 0: get the set of markers
        set<string> mrks_universe;
        for (int ii = 0; ii < num_lgs; ii++) {
            set<string>* tmp_set;
            tmp_set = raw_lgs_[ii]->get_mrks_set();
            mrks_universe.insert(tmp_set->begin(), tmp_set->end());
            delete tmp_set;
        }
        
        map<string, vector<int> > marker_signatures;
        for (set<string>::iterator iter1=mrks_universe.begin();
             iter1 != mrks_universe.end();
             ++iter1) {
            marker_signatures[*iter1] = vector<int>(num_lgs, -1);
        }
        
        // Step 1: Generate the signatures for all the markers
        for (int ii = 0; ii < num_lgs; ii++) {
            // for each linkage group, update its markers signature
            const vector<vector<string> >& mrks_bins = raw_lgs_[ii]->get_markers_in_bins();
            int num_bins = mrks_bins.size();
            for (int jj = 0; jj < num_bins; jj++) {
                for (vector<string>::const_iterator iter2 = mrks_bins[jj].begin();
                     iter2 != mrks_bins[jj].end();
                     ++iter2) {
                     marker_signatures[*iter2][ii] = jj;
                }
            }
        }
        
        // Step 2: Now update every LG to consolidate its markers
        for (int ii = 0; ii < num_lgs; ii++) {
            double lower_bound = raw_lgs_[ii]->get_lower_bound();
            double upper_bound = raw_lgs_[ii]->get_upper_bound();
            int number_of_bins = raw_lgs_[ii]->get_number_of_bins();
            string description = raw_lgs_[ii]->get_description();
            double confidence = raw_lgs_[ii]->get_confidence();
            vector<double> distances = raw_lgs_[ii]->get_distance_btw_adjacent_pairs();
            vector<vector<string> > mrks_bins(number_of_bins);
            int number_of_mrks = 0;
            const vector<vector<string> >& old_mrks_bins = raw_lgs_[ii]->get_markers_in_bins();
            for (int jj = 0; jj < old_mrks_bins.size(); jj++) {
                set<string> crt_bin_mrks;
                crt_bin_mrks.insert(old_mrks_bins[jj].begin(), old_mrks_bins[jj].end());
                while (crt_bin_mrks.size() > 0) {
                    string crt_mrk = *(crt_bin_mrks.begin()); 
                    set<string> crt_sn;
                    crt_sn.insert(crt_mrk);
                    for (set<string>::iterator iter3 = crt_bin_mrks.begin();
                         iter3 != crt_bin_mrks.end();
                         ++iter3) {
                        if (IdenticalSig(marker_signatures[crt_mrk], marker_signatures[*iter3])) {
                            crt_sn.insert(*iter3);
                        }
                    }

                    // delete crt_sn from crt_bin_mrks
                    string crt_sn_str = "";
                    for (set<string>::iterator iter4 = crt_sn.begin(); iter4 != crt_sn.end(); ++iter4) {
                        crt_bin_mrks.erase(*iter4);
                        crt_sn_str = crt_sn_str + *iter4 + ",";
                    }

                    crt_sn_str = crt_sn_str.substr(0, crt_sn_str.length() - 1);
                    mrks_bins[jj].push_back(crt_sn_str);
                    number_of_mrks = number_of_mrks + 1;
                }
            }
            linkage_group_bin* tmp_ptr = new linkage_group_bin();
            tmp_ptr->initialize(number_of_bins,
                                number_of_mrks,
                                lower_bound,
                                upper_bound,
                                description,
                                confidence,
                                distances,
                                mrks_bins);
            delete raw_lgs_[ii];
            raw_lgs_[ii] = tmp_ptr;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    FragClstr* LG_CLUSTER::GenFragClstr(const vector<string>& scc){
        vector<Fragment> frags;
        vector<int> ids;
        for (int ii = 0; ii < raw_lgs_.size(); ii++) {
            Fragment* tmp_frag = raw_lgs_[ii]->GenFrag(scc);
            if (tmp_frag != NULL) {
                if (tmp_frag->get_num_bins() > 1) {
                    frags.push_back(*tmp_frag);
                    ids.push_back(ii);
                }
                delete tmp_frag;
            }
        }
        FragClstr* to_return = new FragClstr();
        to_return->Initialize(frags, ids);
        return to_return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    SCC_Drawer* LG_CLUSTER::GenDrawer(int scc_id, bool direction, vector<string>& scc) {
        assert(scc.size() > 0);
        if (scc.size() == 1) {
            SCC_Drawer* to_return = new Sigleton_Drawer(scc[0], scc_id);
            return to_return;
        } else {
            FragClstr* tmp_clstr = this->GenFragClstr(scc);
            if (direction == false) { // to draw in the opposite direction
                tmp_clstr->FlipOrientation();
            }
            tmp_clstr->solveLP(LP_epsilon);
            // record the markers to be deleted
            vector<pair<string, int> > markers_to_delete = tmp_clstr->return_conflicts();
            accumulate_conflicts(markers_to_delete);
            
            SCC_Drawer* to_return = new Multi_Drawer(*tmp_clstr, scc_id);
            assert(tmp_clstr != NULL);
            delete tmp_clstr;
            return to_return;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLUSTER::gen_consensus_LG(){
        /*
            It is assumed that the orientation for the individual LGs has already been fixed
        */
        int num_lgs = raw_lgs_.size();
        
        // The consensus_lgs_dag needs to be deleted at the end
        linkage_group_DAG* consensus_lgs_dag = NULL;
        
        // Step 0: construct the linkage_group_DAG object
        assert(raw_lgs_.size() > 0);
        consensus_lgs_dag = raw_lgs_[0]->construct_linkage_group_DAG(true);

        for (int mm = 1; mm < num_lgs; mm++) {
            linkage_group_DAG* tmp_lg_dag1 = raw_lgs_[mm]->construct_linkage_group_DAG(true);
            linkage_group_DAG* tmp_lg_dag2 = consensus_lgs_dag;            
            consensus_lgs_dag = combine_lg_DAGs(tmp_lg_dag1, tmp_lg_dag2);
            delete tmp_lg_dag1;
            delete tmp_lg_dag2;
        }
        
        consensus_LG = consensus_lgs_dag->construct_linkage_group_graph_representation();
        assert(consensus_LG != NULL); // crash the program if something gone bad
        // Clean up
        assert(consensus_lgs_dag != NULL);
        delete consensus_lgs_dag;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void LG_CLUSTER::Dump() {
        cout << "# of markers:" << consensus_LG->get_num_of_markers()
             << " #strongly connected components:"<< consensus_LG->get_num_of_scc()
             << " # singletons:" << consensus_LG->get_num_of_singletons()
             << " # redundent edges:" << consensus_LG->get_total_number_of_redundent_edges()
             << endl;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLUSTER::Gen_dot_files(string file_path){
        ofstream output_file(file_path.c_str());
        string to_output = consensus_LG->output_dot_scc(this);
        output_file << to_output;
        output_file.close();
        // set conflicts_computed to true after the first call of this function
        conflicts_computed = true;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void LG_CLUSTER::accumulate_conflicts(const vector<pair<string, int> >& conflicts) {
        assert(conflicts_computed == false);
        for (int ii = 0; ii < conflicts.size(); ii++) {
            string mrkr = conflicts[ii].first;
            int frag_id = conflicts[ii].second;
            assert(frag_id >= 0);
            assert(frag_id < raw_lgs_.size());
            nodes_to_remove[frag_id].push_back(mrkr);
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    void LG_CLUSTER::remove_conflicts() {
        assert(conflicts_computed);
        for (int ii = 0; ii < raw_lgs_.size(); ii++) {
            raw_lgs_[ii]->remove_markers(nodes_to_remove[ii]);
        };
        
        // now re-compute the consensus_LG
        assert(consensus_LG != NULL);
        condense_into_sn_scc();
        condense_into_sn();
        delete consensus_LG;
        consensus_LG = NULL;
        gen_consensus_LG();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    void LG_CLUSTER::Linearize_Dag() {
        consensus_LG->compute_linear_order();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    void LG_CLUSTER::Gen_linear_map(string file_path) {
        ofstream output_file(file_path.c_str());
        string to_output = consensus_LG->output_linear_order();
        output_file << to_output;
        output_file.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    string LG_CLUSTER::Gen_map_chart(){
        return consensus_LG->output_linear_map();
    };

    LG_CLSTR_ENS::LG_CLSTR_ENS() {
        return;
    }
    
    void LG_CLSTR_ENS::Initialize(vector<LG_CLUSTER*>& _lg_clusters) {
        for (int i = 0; i < _lg_clusters.size(); i++){
            assert(_lg_clusters[i] != NULL);
            lg_clusters_.push_back(_lg_clusters[i]);
            _lg_clusters[i] = NULL;
        }
    };
    
    void LG_CLSTR_ENS::Gen_dot_files(){
        for (int i = 0; i < lg_clusters_.size(); i++) {
            char file_name[FILE_NAME_MAX_LENGTH];
            sprintf(file_name,"lg%d.dot",i);
            lg_clusters_[i]->Gen_dot_files(file_name);
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void LG_CLSTR_ENS::Gen_dot_conflict_free() {
        for (int i = 0; i < lg_clusters_.size(); i++) {
            char file_name[FILE_NAME_MAX_LENGTH];
            sprintf(file_name,"lg%d_consensus.dot",i);
            lg_clusters_[i]->Gen_dot_files(file_name);
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void LG_CLSTR_ENS::Gen_linear_graph() {
        for (int i = 0; i < lg_clusters_.size(); i++) {
            char file_name[FILE_NAME_MAX_LENGTH];
            sprintf(file_name,"lg%d_linear.dot",i);
            lg_clusters_[i]->Gen_linear_map(file_name);
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void LG_CLSTR_ENS::Gen_linear_map(){ // gennerate the map in mapchart format
        ofstream output_file("linear_map_chart.txt");
        for (int i = 0; i < lg_clusters_.size(); i++) {
            char buffer[30];
            sprintf(buffer, "group lg%d", i);
            output_file << endl << endl;
            output_file << buffer << endl;
            string lg_map = lg_clusters_[i]->Gen_map_chart();
            output_file << lg_map;
        }
        output_file.close();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
    LG_CLSTR_ENS::~LG_CLSTR_ENS() {
        for (int i = 0; i < lg_clusters_.size(); i++) {
            assert(lg_clusters_[i] != NULL);
            delete lg_clusters_[i];
        }
    };
}
