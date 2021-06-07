/*
 *  shortest_path.cpp
 *  consensusmap
 *
 *  Created by yonghui on 11/9/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "shortest_path.h"
#include <limits>
#include <iostream>
#include <stdio.h>
#include <cmath>
#include <map>
#include <cstdlib>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
using namespace boost;

namespace consensus_map {
    CG::~CG(){
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

    void CG::solve_lp(double epsilon){
        bool no_delete = valid_solution();
        if (no_delete) { // if there is no need to delete nodes
            cout << "CG::solve_lp() no need to delete nodes" << endl;
            return;
        }
        cout << "problem size: " << n_vertex_ << endl;
        double max_weight = 0;
        for (int ii = 0; ii < n_vertex_; ii++) {
            if (g_[ii].weight > max_weight) {
                max_weight = g_[ii].weight;
            }
        }
        // find the initial solution, which is the longest path in terms of the weights assigned
        for (int ii = 0; ii < n_vertex_; ii++) {
            g_[ii].dual_val = 1 + max_weight - g_[ii].weight;
        }
        double dummy1;
        vector<int> initial_path;
        find_shortest_path(dummy1, initial_path);
        for (vector<int>::iterator iter1 = initial_path.begin(); iter1 != initial_path.end(); ++iter1) {
            g_[*iter1].total_flow = 1;
        }
        
        double lambda = compute_lambda();
        double alpha = (4 / (lambda * epsilon)) * log(2 * n_vertex_ / epsilon);
        double delta = epsilon / (4 * alpha);

        double crt_left = 0.0;
        double crt_lambda = lambda;
        int iteration = -1;

        while (true) {
            iteration++;
            compute_dual_vals(alpha);

            double middle = 0; // middle corresponds to \vec{z}^tA\vec{y}
            for (int ii = 0; ii < n_vertex_; ii++) {
                middle = middle + g_[ii].dual_val * g_[ii].total_flow;
            }
            
            double right = 0; // right corresponds to \lambda\vec{z}^t\vec{b}
            for (int ii = 0; ii < n_vertex_; ii++) {
                right = right + crt_lambda * g_[ii].dual_val * g_[ii].weight;
            }
            
            vector<int> opt_path;
            double left; // left corresponds to \vec{z}^tA\vec{y}
            find_shortest_path(left, opt_path);
            crt_left = left;
            
            if (iteration % 1000 == 0) {
                cout << "iteration " << iteration 
                     << "  left:" << left 
                     << "  middle:" << middle 
                     << "  right:" << right
                     << "  sol:" << 1 / crt_lambda 
                     << endl;
            }

            // check if the following condition is met or not
            // \vec{z}^tA\vec{y} - C(\vec{z}) \le \epsilon(\vec{z}^tA\vec{y} + \lambda\vec{z}^t\vec{b})
            if ((middle - left) < epsilon * (middle + right) ) {
                cout << "converged!" << endl;
                cout << "iteration " << iteration 
                     << "  left:" << left 
                     << "  middle:" << middle 
                     << "  right:" << right
                     << "  sol:" << 1 / crt_lambda 
                     << endl;
                break;
            } else {
                update_solution(delta, opt_path);
                crt_lambda = compute_lambda();
                if (crt_lambda < 0.5 * lambda) {
                    lambda = crt_lambda;
                    alpha = 4 / (lambda * epsilon) * log(2 * n_vertex_ / epsilon);
                    delta = epsilon / (4 * alpha);
                }
            }
        }
        
        // compute the normalized dual solution
        for (int ii = 0; ii < n_vertex_; ii++) {
            g_[ii].nor_dual_val = g_[ii].dual_val / crt_left;
            if (g_[ii].nor_dual_val > 1) {
                g_[ii].nor_dual_val = 1;            
            }
        }
        // as the last step, round the fractional solution to integral solutions
        round_lp();    
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CG::dump(){
        for (int ii = 0; ii < n_vertex_; ii++) {
            printf("weight: %.3f, dual_val: %.3f, total_flow: %.3f, nor_dual_val:%.3f, n_edges:%d, to_delete:%d", 
                   g_[ii].weight,
                   g_[ii].dual_val,
                   g_[ii].total_flow,
                   g_[ii].nor_dual_val,
                   g_[ii].num_edges,
                   g_[ii].to_delete);
            for (int jj = 0; jj < g_[ii].num_edges; jj++) {
                cout << ";" << g_[ii].edges[jj].to_id << " " << g_[ii].edges[jj].special;
            }
            cout << endl;
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    V_Satellite* CG::construct_sub_graph(const vector<int>& vertices) const{
        map<int,int> old_new;
        int sub_size = vertices.size();
        for (int ii = 0; ii < sub_size; ii++) {
            old_new[vertices[ii]] = ii;
        }
        V_Satellite* sub_g = new V_Satellite[sub_size];
        vector<vector<E_Satellite> > edges;
        edges.resize(sub_size);
        for (int ii = 0; ii < sub_size; ii++){
            int node_id = vertices[ii];
            for (int jj = 0; jj < g_[node_id].num_edges; jj++) {
                int to_id = g_[node_id].edges[jj].to_id;
                if (old_new.find(to_id) != old_new.end()) {
                    int new_to_id = old_new[to_id];
                    bool is_special = g_[node_id].edges[jj].special;
                    E_Satellite new_edge = {new_to_id, is_special};
                    edges[ii].push_back(new_edge);
                }
            }
        }
        for (int ii = 0; ii < sub_size; ii++) {
            sub_g[ii] = g_[vertices[ii]];
            sub_g[ii].dual_val = 0.0;
            sub_g[ii].total_flow = 0.0;
            sub_g[ii].nor_dual_val = 0.0;
            sub_g[ii].to_delete = false;
            sub_g[ii].num_edges = edges[ii].size();
            sub_g[ii].edges = new E_Satellite[edges[ii].size()];
            for (int jj = 0; jj < edges[ii].size(); jj++) {
                sub_g[ii].edges[jj] = edges[ii][jj];
            }
        }
        
        return sub_g;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CG::compute_dual_vals(double alpha) {
        // compute on the log-scale first to avoid overflow
        for (int ii = 0; ii < n_vertex_; ii++) {
            double dual_val = log(1 / g_[ii].weight) + (alpha * g_[ii].total_flow / g_[ii].weight);
            g_[ii].dual_val = dual_val;
        }
        double max_dual = 0;
        for (int ii = 0; ii < n_vertex_; ii++) {
            if (g_[ii].dual_val > max_dual) {
                max_dual = g_[ii].dual_val;
            }
        }
        for (int ii = 0; ii < n_vertex_; ii++) {
            g_[ii].dual_val = g_[ii].dual_val - max_dual + 50;
            g_[ii].dual_val = exp(g_[ii].dual_val);
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    double CG::compute_lambda() {
        double max_lambda = 0;
        for (int ii = 0; ii < n_vertex_; ii++) {
            double lambda = g_[ii].total_flow / g_[ii].weight;
            if (lambda > max_lambda) {
                max_lambda = lambda;
            }
        }
        return max_lambda;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void CG::update_solution(double delta, vector<int> & path) {
        for (int ii = 0; ii < n_vertex_; ii++) {
            g_[ii].total_flow = (1 - delta) * g_[ii].total_flow;
        }
        for (vector<int>::iterator iter1 = path.begin(); iter1 != path.end(); ++iter1) {
            g_[*iter1].total_flow = g_[*iter1].total_flow + delta;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    void CG::find_shortest_path(double& crt_cost, vector<int>& opt_path){
        vector<int> crt_shortest_path;
        double cost = numeric_limits<double>::max();
        
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
            }
            
            double min_dist1 = numeric_limits<double>::max();
            int closest_dest = -1;
            for (vector<int>::iterator sink = sinks.begin();
                 sink != sinks.end();
                 ++sink) {
                if (node_status[*sink].distance <= min_dist1) {
                    min_dist1 = node_status[*sink].distance;
                    closest_dest = *sink;
                }
            }
            assert(closest_dest != -1);
            if (min_dist1 < cost) {
                cost = min_dist1;
                crt_shortest_path.clear();
                crt_shortest_path.push_back(closest_dest);
                int parent = node_status[closest_dest].parent;
                while (parent > -1) {
                    crt_shortest_path.push_back(parent);
                    parent = node_status[parent].parent;
                }
            }
                        
        } // end of for (dijastra to find shortest path from each node)
        
        crt_cost = cost;
        assert(cost < numeric_limits<double>::max());
        delete[] node_status;
        opt_path = crt_shortest_path;
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CG::initialize_lp_round() {
        for (int ii = 0; ii < n_vertex_; ii++) {
            g_[ii].to_delete = false;
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void CG::round_lp_randomized(vector<int> nodes){ 
        // step0 decide if the sub-problem is solved already

        map<int, int> ori_2_new_map;
        for (int ii = 0; ii < nodes.size(); ii++) {
            ori_2_new_map[nodes[ii]] = ii;
        }
        
        double total_probs = 0.0; // care should be taken when dealing with nodes that have been deleted already
        vector<double> probs;
        probs.resize(nodes.size());
        for (int ii = 0; ii < nodes.size(); ii++) {
            if (g_[nodes[ii]].to_delete) {
                probs[ii] = 0.0;
            } else {
                probs[ii] = g_[nodes[ii]].nor_dual_val;
                total_probs = total_probs + probs[ii];
            }
        }
        
        bool exist_none_special_edges = exist_no_special_edge(nodes);
        
        if (not exist_none_special_edges) {
            return;
        }

        // step1: sample one node to delete
        // random_double is a number between 0 and total_probs        
        double random_double = ((rand() + 0.5) / (double(RAND_MAX) + 1)) * total_probs;
        cout << "dbg print:random_double" << random_double / total_probs << endl;
        double accumulated_probs = 0.0;
        int delete_index;
        for (int ii = 0; ii < nodes.size(); ii++) {
            delete_index = ii;
            if (probs[ii] == 0.0) continue;
            accumulated_probs = accumulated_probs + probs[ii];
            if (accumulated_probs > random_double) break;
        }
        
        // step2: delete the vertex
        g_[nodes[delete_index]].to_delete = true;

        
        // step 3: construct a graph. Find out the strongly connected components
        typedef boost::adjacency_list<vecS,vecS,directedS> Graph;
        Graph G(nodes.size());
        for (int ii = 0; ii < nodes.size(); ii++) {
            if (not g_[nodes[ii]].to_delete) {
                for (int jj = 0; jj < g_[nodes[ii]].num_edges; jj++) {
                    int to_id = g_[nodes[ii]].edges[jj].to_id;
                    if ((ori_2_new_map.find(to_id) != ori_2_new_map.end()) and 
                        (not g_[to_id].to_delete)) {
                        add_edge(ii, ori_2_new_map[to_id], G);                    
                    }
                }
            }
        }        
        vector<int> components(nodes.size(),-1);
        int num = strong_components(G,&components[0]);

        // step4: recurse
        vector<vector<int> > sccs;
        sccs.resize(num);
        for (int ii = 0; ii < nodes.size(); ii++){
            assert(components[ii] < num);
            assert(components[ii] >= 0);
            sccs[components[ii]].push_back(nodes[ii]);
        }
        for (int ii = 0; ii < num; ii++) {
            round_lp_randomized(sccs[ii]);
        }
        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void CG::round_lp_greedy(vector<int> nodes){ 
        // step0 decide if the sub-problem is solved already

        map<int, int> ori_2_new_map;
        for (int ii = 0; ii < nodes.size(); ii++) {
            ori_2_new_map[nodes[ii]] = ii;
        }
        
        double total_probs = 0.0; // care should be taken when dealing with nodes that have been deleted already
        vector<double> probs;
        probs.resize(nodes.size());
        for (int ii = 0; ii < nodes.size(); ii++) {
            if (g_[nodes[ii]].to_delete) {
                probs[ii] = 0.0;
            } else {
                probs[ii] = g_[nodes[ii]].nor_dual_val;
                total_probs = total_probs + probs[ii];
            }
        }
        
        bool exist_none_special_edges = exist_no_special_edge(nodes);
        
        if (not exist_none_special_edges) {
            return;
        }

        // step1: sample one node to delete
        // use a greedy approach to decide which node to delete
        int delete_index = -1;
        double max_prob = 0.0;
        for (int ii = 0; ii < nodes.size(); ii++) {
            if (probs[ii] > max_prob) {
                delete_index = ii;
                max_prob = probs[ii];
            }
        }
        assert(delete_index != -1);
        
        // step2: delete the vertex
        g_[nodes[delete_index]].to_delete = true;

        
        // step 3: construct a graph. Find out the strongly connected components
        typedef boost::adjacency_list<vecS,vecS,directedS> Graph;
        Graph G(nodes.size());
        for (int ii = 0; ii < nodes.size(); ii++) {
            if (not g_[nodes[ii]].to_delete) {
                for (int jj = 0; jj < g_[nodes[ii]].num_edges; jj++) {
                    int to_id = g_[nodes[ii]].edges[jj].to_id;
                    if ((ori_2_new_map.find(to_id) != ori_2_new_map.end()) and 
                        (not g_[to_id].to_delete)) {
                        add_edge(ii, ori_2_new_map[to_id], G);                    
                    }
                }
            }
        }        
        vector<int> components(nodes.size(),-1);
        int num = strong_components(G,&components[0]);

        // step4: recurse
        vector<vector<int> > sccs;
        sccs.resize(num);
        for (int ii = 0; ii < nodes.size(); ii++){
            assert(components[ii] < num);
            assert(components[ii] >= 0);
            sccs[components[ii]].push_back(nodes[ii]);
        }
        for (int ii = 0; ii < num; ii++) {
            round_lp_greedy(sccs[ii]);
        }
        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    bool CG::exist_no_special_edge(const vector<int>& nodes){
        map<int, int> ori_2_new_map;
        for (int ii = 0; ii < nodes.size(); ii++) {
            ori_2_new_map[nodes[ii]] = ii;
        }
        bool exist_none_special_edges = false;
        for (int ii = 0; ii < nodes.size(); ii++) {
            if (not g_[nodes[ii]].to_delete) {
                for (int jj = 0; jj < g_[nodes[ii]].num_edges; jj++) {
                    int to_id = g_[nodes[ii]].edges[jj].to_id;
                    if ((ori_2_new_map.find(to_id) != ori_2_new_map.end()) and 
                        (not g_[to_id].to_delete)) {
                        if (not g_[nodes[ii]].edges[jj].special) {
                            exist_none_special_edges = true;
                        }
                    }
                }
            }
        }
        return exist_none_special_edges;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool CG::valid_solution(){
        typedef boost::adjacency_list<vecS, vecS, directedS> Graph;
        Graph G(n_vertex_);
        for (int ii = 0; ii < n_vertex_; ii++) {
            if (not g_[ii].to_delete) {
                for (int jj = 0; jj < g_[ii].num_edges; jj++) {
                    int to_id = g_[ii].edges[jj].to_id;
                    if (not g_[to_id].to_delete) {
                        add_edge(ii, to_id, G);                    
                    }
                }
            }
        }        
        vector<int> components(n_vertex_, -1);
        int num = strong_components(G, &components[0]);

        // check each connected component being a simple component
        vector<vector<int> > sccs;
        sccs.resize(num);
        for (int ii = 0; ii < n_vertex_; ii++){
            assert(components[ii] < num);
            assert(components[ii] >= 0);
            sccs[components[ii]].push_back(ii);
        }
        for (int ii = 0; ii < num; ii++) {
            bool none_special_edges = exist_no_special_edge(sccs[ii]);
            if (none_special_edges) {
                return false;
            }
        } 
        return true;       
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    void CG::minimize_sol(){
        bool exist_redundent = true;
        int count = 0;
        while (exist_redundent) {
            exist_redundent = false;
            for (int ii = 0; ii < n_vertex_; ii++) {
                if (g_[ii].to_delete) {
                    // test to put back the node
                    g_[ii].to_delete = false;
                    bool valid = valid_solution();
                    if (valid) {
                        exist_redundent = true;
                        count++;
                    } else {
                        g_[ii].to_delete = true;
                    }
                }
            }
        }
        cout << "removed " << count << " redundent nodes" << endl;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    void CG::round_lp() {
        vector<int> entire_set;
        for (int ii = 0; ii < n_vertex_; ii++) {
            entire_set.push_back(ii);
        }
        initialize_lp_round();
        round_lp_greedy(entire_set);
        minimize_sol();
        vector<int> greedy_sol;
        for (int ii = 0; ii < n_vertex_; ii++) {
            if (g_[ii].to_delete) {
                greedy_sol.push_back(ii);
            }
        }
        initialize_lp_round();
        round_lp_randomized(entire_set);
        minimize_sol();
        vector<int> randomized_sol;
        for (int ii = 0; ii < n_vertex_; ii++) {
            if (g_[ii].to_delete) {
                randomized_sol.push_back(ii);
            }
        }
        
        cout << "randomized solution size:" << randomized_sol.size() << endl;
        cout << "greedy solution size:" << greedy_sol.size() << endl;
        initialize_lp_round();
        if (greedy_sol.size() < randomized_sol.size()) {
            for (int ii = 0; ii < greedy_sol.size(); ii++) {
                g_[greedy_sol[ii]].to_delete = true;
            }
        } else {
            for (int ii = 0; ii < randomized_sol.size(); ii++) {
                g_[randomized_sol[ii]].to_delete = true;
            }
        }
        
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
