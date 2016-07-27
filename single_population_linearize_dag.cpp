/*
 *  single_population_linearize_dag.cpp
 *  consensusmap
 *
 *  Created by yonghui on 11/30/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#include <cassert>
#include "single_population.h"
#include "boost/graph/topological_sort.hpp"
#include <queue>

using namespace std;

namespace consensus_map{

    void print_vector(const vector<int> & v_in ) {
        for (int ii = 0; ii < v_in.size(); ii++) {
            cout << v_in[ii] << '\t';
        }
        cout << endl;
    }
    void linkage_group_graph_representation::initialize_linearization() {
        // make sure the graph is a dag
        assert(num_scc == number_of_markers);
        linear_order.resize(number_of_markers);
        pair_wise_distances.resize(number_of_markers);
        for (int ii = 0; ii < number_of_markers; ii++) {
            pair_wise_distances[ii].resize(number_of_markers);
            for (int jj = 0; jj < number_of_markers; jj++) {
                pair_wise_distances[ii][jj] = -1;
            }
        }
        for (int ii = 0; ii < number_of_markers; ii++) {
            pair_wise_distances[ii][ii] = 0;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   
    void linkage_group_graph_representation::compute_pair_wise_distances() {
        // do topological sort first
        vector<int> in_degree(number_of_markers, 0);
        vector<bool> visited(number_of_markers, false);
        vector<int> tp_order;

        for (int ii = 0; ii < number_of_markers; ii++) {
            for (int jj = 0; jj < adjacency_list[ii].size(); jj++) {
                int to_id = adjacency_list[ii][jj].to_id;
                in_degree[to_id] = in_degree[to_id] + 1;
            }
        };
        
        priority_queue<TP_Sort_NS, vector<TP_Sort_NS>, compare_NS> pq;
        for (int ii = 0; ii <number_of_markers; ii++) {
            TP_Sort_NS ns = {ii, in_degree[ii]};
            pq.push(ns);
        }
        
        while (not pq.empty()) {
            TP_Sort_NS ns = pq.top();
            pq.pop();
            int cur_id = ns.node_id;
            if (not visited[cur_id]) {
                assert(in_degree[cur_id] == 0);
                tp_order.push_back(cur_id);
                visited[cur_id] = true;
                for (int jj = 0; jj < adjacency_list[cur_id].size(); jj++) {
                    int to_id = adjacency_list[cur_id][jj].to_id;
                    in_degree[to_id] = in_degree[to_id] - 1;
                    TP_Sort_NS ns = {to_id, in_degree[to_id]};
                    pq.push(ns);
                }
            }
        }
        assert(tp_order.size() == number_of_markers);
        tp_linear_order = tp_order;
        // Now compute the pairwise distance
        // The order of the computation is important, start from the last element to the first 
        // element in tp order
        for (int ii = number_of_markers - 1; ii >= 0; ii--) {
            int from_id = tp_order[ii];
            for (int jj = 0; jj < number_of_markers; jj++) {
                int count = 0;
                double total = 0;
                for (int kk = 0; kk < adjacency_list[from_id].size(); kk++) {
                    int to_id = adjacency_list[from_id][kk].to_id;
                    if (pair_wise_distances[to_id][jj] >= 0) {
                        total = total + pair_wise_distances[to_id][jj];
                        total = total + adjacency_list[from_id][kk].weight;
                        count = count + 1;
                    }
                }
                if (count > 0) {
                    pair_wise_distances[from_id][jj] = total / count;
                }
            }
        }
        
        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void linkage_group_graph_representation::impute_missing_distance() {
        for (int ii = 0; ii < number_of_markers; ii++) {
            for (int jj = ii + 1; jj < number_of_markers; jj++) {
                if ((pair_wise_distances[ii][jj] < 0) and (pair_wise_distances[jj][ii] < 0)) {
                    // find the set of common parents and common children
                    vector<int> common_parents;
                    vector<int> common_children;
                    for (int kk = 0; kk < number_of_markers; kk++) {
                        if ((pair_wise_distances[kk][ii] > 0) and 
                            (pair_wise_distances[kk][jj] > 0) and 
                            (kk != ii) and 
                            (kk != jj)){
                            common_parents.push_back(kk);
                        }
                        if ((pair_wise_distances[ii][kk] > 0) and 
                            (pair_wise_distances[jj][kk] > 0) and 
                            (kk != ii) and 
                            (kk != jj)){
                            common_children.push_back(kk);
                        }
                    }
                    assert((common_parents.size() > 0)or (common_children.size() > 0));
                    if ((common_parents.size() > 0) and (common_children.size() > 0)) {
                        double distance = 0;
                        for (int mm = 0; mm < common_parents.size(); mm++) {
                            for (int nn = 0; nn < common_children.size(); nn++) {
                                int parent_id = common_parents[mm];
                                int child_id = common_children[nn];
                                double ii_frac = pair_wise_distances[parent_id][ii] 
                                    / (pair_wise_distances[parent_id][ii] + pair_wise_distances[ii][child_id]);
                                double jj_frac = pair_wise_distances[parent_id][jj] 
                                    / (pair_wise_distances[parent_id][jj] + pair_wise_distances[jj][child_id]);
                                distance = distance + (jj_frac - ii_frac) * pair_wise_distances[parent_id][child_id];
                            }
                        }
                        distance = distance / (common_parents.size() * common_children.size());
                        pair_wise_distances[ii][jj] = distance;
                        pair_wise_distances[jj][ii] = -distance;
                    } else if (common_parents.size() > 0) {
                        double distance = 0;
                        for (int mm = 0; mm < common_parents.size(); mm++) {
                            int parent_id = common_parents[mm];
                            distance = distance + pair_wise_distances[parent_id][jj] - pair_wise_distances[parent_id][ii];
                        }
                        distance = distance / common_parents.size();
                        pair_wise_distances[ii][jj] = distance;
                        pair_wise_distances[jj][ii] = -distance;
                    } else if (common_children.size() > 0) {
                        double distance = 0;
                        for (int nn = 0; nn < common_children.size(); nn++) {
                            int child_id = common_children[nn];
                            distance = distance + pair_wise_distances[ii][child_id] - pair_wise_distances[jj][child_id];
                        }
                        distance = distance / common_children.size();
                        pair_wise_distances[ii][jj] = distance;
                        pair_wise_distances[jj][ii] = -distance;
                    }
                }
            }
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int linkage_group_graph_representation::select_node(const vector<int>& top_nodes){
        if (top_nodes.size() == 1) {
            return top_nodes[0];
        }
        double min_weight = numeric_limits<double>::max();
        int min_index = -1;
        // find a node with the lowest in-degree 
        for (int ii = 0; ii < top_nodes.size(); ii++) {
            double in_weights = 0;
            int from_id = top_nodes[ii];
            for (int jj = 0; jj < top_nodes.size(); jj++) {
                if (ii != jj) {
                    int to_id = top_nodes[jj];
                    if (pair_wise_distances[to_id][from_id] > 0) {
                        in_weights = in_weights + pair_wise_distances[to_id][from_id];
                    }
                }
            }
            if (in_weights < min_weight) {
                min_weight = in_weights;
                min_index = ii;
            }
        }
        assert(min_index >= 0);
        if ( min_weight > 0) {
            cout << "min weight:" << min_weight << endl;
        }
        return top_nodes[min_index];
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void linkage_group_graph_representation::linearize_graph() {
        // do topological sort first
        vector<int> in_degree(number_of_markers, 0);
        vector<bool> visited(number_of_markers, false);
        vector<int> tp_order;

        for (int ii = 0; ii < number_of_markers; ii++) {
            for (int jj = 0; jj < adjacency_list[ii].size(); jj++) {
                int to_id = adjacency_list[ii][jj].to_id;
                in_degree[to_id] = in_degree[to_id] + 1;
            }
        };
        
        for (int ii = 0; ii < number_of_markers; ii++) {
            vector<int> top_nodes;
            for (int jj = 0; jj < number_of_markers; jj++) {
                if ((in_degree[jj] == 0) and (not visited[jj])) {
                    top_nodes.push_back(jj);
                }
            }
            assert(top_nodes.size() >= 1);
            int next_node = select_node(top_nodes);
            tp_order.push_back(next_node);
            visited[next_node] = true;
            for (int kk = 0; kk < adjacency_list[next_node].size(); kk++) {
                int to_id = adjacency_list[next_node][kk].to_id;
                in_degree[to_id] = in_degree[to_id] - 1;
            }
        }
        linear_order = tp_order;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    string linkage_group_graph_representation::output_linear_order() {
        // To draw in the forward direction
        
        ostringstream ss_out(ios_base::out);
        ss_out << "digraph G {" << endl;
        
        // draw the edges
        for (int ii = 0; ii < number_of_markers-1; ii++) {
            double weight = pair_wise_distances[linear_order[ii]][linear_order[ii+1]];
            char buffer[100];
            sprintf(buffer, "%.2f",weight);
            ss_out << "\"" << marker_ids[linear_order[ii]] << "\"" 
                    << " -> " 
                    << "\"" << marker_ids[linear_order[ii+1]] << "\"" 
                    << "[label =" << buffer << "]" << ';' << endl;
        }        
        
        ss_out << "}" << endl << endl;

        return ss_out.str();
    } 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    string linkage_group_graph_representation::output_linear_map(){ // output the final linearized map in mapchart format
        ostringstream ss_out(ios_base::out);
        assert(marker_ids.size() > 0);
        ss_out << marker_ids[linear_order[0]] << '\t' << 0.00 << endl;
        double distance = 0;
        for (int ii = 1; ii < number_of_markers; ii++) {
            distance = distance + pair_wise_distances[linear_order[ii - 1]][linear_order[ii]];
            char buffer[40];
            sprintf(buffer, "%.2f", distance);
            ss_out << marker_ids[linear_order[ii]] << '\t' << buffer << endl;
        }
        return ss_out.str();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    void linkage_group_graph_representation::compute_linear_order(){
        initialize_linearization();
        compute_pair_wise_distances();
        impute_missing_distance();
        linearize_graph();        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
