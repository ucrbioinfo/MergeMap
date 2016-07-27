/*
 *  fragment.cpp
 *  consensusmap
 *
 *  Created by yonghui on 10/15/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "fragment.h"
#include "heuristic_solver.h"
#include <iostream>
#include <set>
#include "constants.h"
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>


namespace consensus_map{
    void Fragment::FlipOrientation() {
        int num_bins = markers_in_bins_.size();
        vector<vector<string> > tmp_bins = markers_in_bins_;
        for (int ii = 0; ii < num_bins; ii++) {
            markers_in_bins_[ii] = tmp_bins[num_bins - ii - 1];
        } 
        
        int num_dists = distance_btw_adjacent_pairs_.size();
        vector<double> tmp_dists = distance_btw_adjacent_pairs_;
        for (int ii = 0; ii < num_dists; ii++) {
            distance_btw_adjacent_pairs_[ii] = tmp_dists[num_dists - ii -1];
        }
    }

    void Fragment::dump() {
        for (int ii = 0; ii < markers_in_bins_.size(); ii++) {
            for (int jj = 0; jj < markers_in_bins_[ii].size(); jj++) {
                cout << markers_in_bins_[ii][jj] << "|";
            }
            cout << endl;
        }
        cout << "+++++++++++" << endl;
        for (int ii = 0; ii < prob_del_.size(); ii++) {
            for (int jj = 0; jj < prob_del_[ii].size(); jj++) {
                cout << prob_del_[ii][jj] << '\t';
            }
            cout << endl;
        }

        for (int ii = 0; ii < dels_.size(); ii++) {
            for (int jj = 0; jj < dels_[ii].size(); jj++) {
                cout << dels_[ii][jj] << '\t';
            }
            cout << endl;
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

    vector<pair<string,int> > FragClstr::return_conflicts() {
        vector<pair<string,int> > to_return;
        for (int ii = 0; ii < frags_.size(); ii++) {
            vector<string> to_remove = frags_[ii].return_deleted();
            int frag_id = frag_ids_[ii];
            for (int jj = 0; jj < to_remove.size(); jj++) {
                to_return.push_back(make_pair(to_remove[jj], frag_id));
            }
        }
        return to_return;
    };    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void FragClstr::solveLP(double epsilon) {
        // modified by yonghui on Jan 19 2008
        // solveLP is a misnomer, since this algorithm will try to delete markers based on heuristic algorithms also
        // The main purpose here is to break a graph into a set of strongly connected component 
        // If a strongly connected component contains more than kLargestSCCLP number of markers,
        // a heuristic algorithm is applied first. 
        map<int, node_id_strct> node_id_map;
        vector<vector<int> > sccs;
        int num_scc;
        CG* gg;
        while (true) {
            gg = GenGraphLag(node_id_map);
            int number_of_markers = gg->get_size();
            const V_Satellite* p_vertex = gg->get_graph();
            // break the graph into a set of strongly connected components
            typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> Graph;
            Graph G(number_of_markers);
            
            for (int ii = 0; ii < number_of_markers; ii++) {
                int num_edges_ii = p_vertex[ii].num_edges;
                const E_Satellite* edges_ii = p_vertex[ii].edges;
                for (int jj = 0; jj < num_edges_ii; jj++) {
                    add_edge(ii, edges_ii[jj].to_id, G);
                }
            }
            vector<int> components(number_of_markers, -1);
            num_scc = strong_components(G, &components[0]);
            sccs.clear();
            sccs.resize(num_scc);
            for (int ii = 0; ii < number_of_markers; ii++)
            {
                if (components[ii] >= num_scc){
                    cout << "ERROR! strong component id is unexpected" << endl;
                    assert(false); // crash the program on error
                }
                else {
                    sccs[components[ii]].push_back(ii);
                }
            }
            
            // now for every SCC, check if its size is less than the threshold. If not, call H_Solver to remove one node
            bool one_more_iteration = false;
            for (int ii = 0; ii < num_scc; ii++) {
                if (sccs[ii].size() > kLargestSCCLP) {
                    cout << "the sub-problem size:" << sccs[ii].size() << endl;
                    V_Satellite* sub_g = gg->construct_sub_graph(sccs[ii]);
                    H_Solver sub_problem;
                    sub_problem.Initialize(sub_g, sccs[ii].size());
                    int remove_sub_id = sub_problem.Remove_one_marker();
                    if (remove_sub_id >= 0) {
                        one_more_iteration = true;
                        int remove_id = sccs[ii][remove_sub_id];
                        // now remove the node from the corresponding fragment
                        node_id_strct remove_node_id = node_id_map[remove_id];
                        h_delete(remove_node_id.frag_id, remove_node_id.level_id, remove_node_id.mid);
                    } else {
                        cout << "nothing to remove" << endl;
                    }
                }
            }
            if (not one_more_iteration){
                break;
            } else {
                cout << "one more iteration" << endl;
            }
        }

        // initialize the probs and dels vector
        vector<vector<vector<double> > > probs;
        vector<vector<vector<bool> > > dels;
        probs.resize(frags_.size());
        dels.resize(frags_.size());
        for (int ii = 0; ii < frags_.size(); ii++) {
            const vector<vector<string> > & markers_bins = frags_[ii].get_bins();
            probs[ii].resize(markers_bins.size());
            dels[ii].resize(markers_bins.size());
            for (int jj = 0; jj < markers_bins.size(); jj++){
                probs[ii][jj] = vector<double>(markers_bins[jj].size(), 0);
                dels[ii][jj] = vector<bool>(markers_bins[jj].size(), false);
            }
        }
        
        for (int ii = 0; ii < num_scc; ii++) {
            V_Satellite* sub_g = gg->construct_sub_graph(sccs[ii]);
            CG sub_problem;
            sub_problem.Initialize(sub_g, sccs[ii].size());
            sub_problem.solve_lp(epsilon);
            const V_Satellite* nodes = sub_problem.get_graph();
            for (int jj = 0; jj < sccs[ii].size(); jj++) {
                node_id_strct id = node_id_map[sccs[ii][jj]];
                double prob_jj = nodes[jj].nor_dual_val;
                bool del_jj = nodes[jj].to_delete;
                probs[id.frag_id][id.level_id][id.mid] = prob_jj;
                dels[id.frag_id][id.level_id][id.mid] = del_jj;
            }
        }

        
        for (int ii = 0; ii < frags_.size(); ii++) {
            frags_[ii].set_prob(probs[ii]);
            frags_[ii].set_del(dels[ii]);
        }
        
        
        delete gg;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void FragClstr::dump() {
        for (int ii = 0; ii < frags_.size(); ii++) {
            frags_[ii].dump();
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    CG* FragClstr::GenGraphLag(map<int, node_id_strct>& node_id_map){
        // First get the list of non-singleton markers
        
        // edited by yonghui on Jan 19
        // consider the case where some of the markers are deleted by our heuristic algorithm 
        
        int number_of_h_del = 0;
        map<string, vector<node_id_strct> > marker_presence;
        for (int ii = 0; ii < frags_.size(); ii++){
            const vector<vector<string> > & markers_bins = frags_[ii].get_bins();
            const vector<vector<bool> > & h_del = frags_[ii].get_heuristic_del();
            for (int jj = 0; jj < markers_bins.size(); jj++) {
                for (int kk = 0; kk < markers_bins[jj].size(); kk++) {
                    string crt_mrk = markers_bins[jj][kk];
                    bool crt_del = h_del[jj][kk];
                    if (not crt_del) { // if the marker is not deleted
                        node_id_strct crt_node_id = {ii,jj,kk};
                        if (marker_presence.find(crt_mrk) == marker_presence.end()) {
                            marker_presence[crt_mrk] = vector<node_id_strct>();
                        }
                        marker_presence[crt_mrk].push_back(crt_node_id);
                    } else {
                        number_of_h_del++;
                    }
                }
            }        
        }
        
        set<node_id_strct, ID_cmp> singletons;
        for (map<string, vector<node_id_strct> >::iterator crt_marker = marker_presence.begin();
             crt_marker != marker_presence.end();
             ++crt_marker) {
            if ((crt_marker->second).size() == 1) {
                singletons.insert((crt_marker->second)[0]);
            }
        }        
        
        
        int n_nodes = total_num_mrks() - singletons.size() - number_of_h_del; // this is suspicious
        V_Satellite* vertices = new V_Satellite[n_nodes];
        map<node_id_strct, int, ID_cmp> id_map_rev; // not all node_id's are present
        vector<vector<pair<int,bool> > > adjacent_list;
        adjacent_list.resize(n_nodes);
        // populate the two id maps
        // populate the weight for all the nodes
        int crt_counter = -1;
        for (int ii = 0; ii < frags_.size(); ii++) {
            double confidence = frags_[ii].get_confidence();
            const vector<vector<string> > & markers_bins = frags_[ii].get_bins();
            const vector<vector<bool> > & h_del = frags_[ii].get_heuristic_del();
            for (int jj = 0; jj < markers_bins.size(); jj++) {
            	for (int kk = 0; kk < markers_bins[jj].size(); kk++) {
            	    node_id_strct crt_node_id = {ii,jj,kk};
                    double crt_del = h_del[jj][kk];
            	    if ((not crt_del) and (singletons.find(crt_node_id) == singletons.end())) {
                        crt_counter++;
            	        node_id_map[crt_counter] = crt_node_id;
                        id_map_rev[crt_node_id] = crt_counter;
                        vertices[crt_counter].weight = confidence;
                        vertices[crt_counter].dual_val = 0;
                        vertices[crt_counter].nor_dual_val = 0;
                        vertices[crt_counter].total_flow = 0;
                        vertices[crt_counter].to_delete = false;
                        char buffer[100];
                        sprintf(buffer, "_%d", ii);
                        vertices[crt_counter].name = markers_bins[jj][kk]+buffer;
            	    }
            	}
            }
        }
        
        assert(crt_counter == n_nodes -1);
        
        // construct edges
        for (int ii = 0; ii < frags_.size(); ii++) {
            const vector<vector<string> > & markers_bins = frags_[ii].get_bins();        
            const vector<vector<bool> > & h_del = frags_[ii].get_heuristic_del();
            for (int jj1 = 0; jj1 < markers_bins.size(); jj1++) {
            	for (int kk1 = 0; kk1 < markers_bins[jj1].size(); kk1++) {
            	    node_id_strct crt_node_id1 = {ii,jj1,kk1};
                    bool crt_del1 = h_del[jj1][kk1];
                    for (int jj2 = jj1 + 1; jj2 < markers_bins.size(); jj2++){
                        for (int kk2 = 0; kk2 < markers_bins[jj2].size(); kk2++) {
                            node_id_strct crt_node_id2 = {ii,jj2,kk2};
                            bool crt_del2 = h_del[jj2][kk2];
                            if ((singletons.find(crt_node_id1) == singletons.end()) and 
                                (singletons.find(crt_node_id2) == singletons.end()) and 
                                (not crt_del1) and 
                                (not crt_del2)) {
                                int id1 = id_map_rev[crt_node_id1];
                                int id2 = id_map_rev[crt_node_id2];
                                adjacent_list[id1].push_back(make_pair(id2, false));
                            }
                        }
                    }
            	}
            }
        }
        
        // construct special edges
        for (map<string, vector<node_id_strct> >::iterator crt_marker = marker_presence.begin();
             crt_marker != marker_presence.end();
             ++crt_marker) {
            // cout << crt_marker->first << (crt_marker->second).size() << endl; 
            if ((crt_marker->second).size() > 1) {
                for (vector<node_id_strct>::iterator node1 = (crt_marker->second).begin();
                     node1 != (crt_marker->second).end();
                     ++node1) {
                    int id1 = id_map_rev[*node1];
                    for (vector<node_id_strct>::iterator node2 = (crt_marker->second).begin();
                         node2 != (crt_marker->second).end();
                        ++node2) {
                        int id2 = id_map_rev[*node2];
                        if (id1 != id2){
                           adjacent_list[id1].push_back(make_pair(id2, true));
                        }
                    }
                }
            }
        }
        
        for (int ii = 0; ii < n_nodes; ii++) {
            int n_edges = adjacent_list[ii].size();
            vertices[ii].num_edges = n_edges;
            vertices[ii].edges = new E_Satellite[n_edges];
            for (int jj = 0; jj < n_edges; jj++) {
                vertices[ii].edges[jj].to_id = adjacent_list[ii][jj].first;
                vertices[ii].edges[jj].special = adjacent_list[ii][jj].second;                
            }
        }
        
        CG* to_return = new CG();
        to_return->Initialize(vertices, n_nodes);        
        return to_return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} // end of the namespace
