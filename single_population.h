/*
 *  single_population.h
 *  consensusmap
 *
 *  Created by yonghui on 5/12/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef SINGLE_POPULATION_HEADER
#define SINGLE_POPULATION_HEADER

#include <vector>
#include <set>
#include <string>
#include <utility>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <sstream>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include "constants.h"
#include "fragment.h"
#include <map>
#include <utility>
#include "ConsensusLG.h"
using namespace std;
using namespace boost;

/*all the classes are defined in the consensus_map namespace*/
namespace consensus_map
{

    // forward declaration
    class LG_CLUSTER;
    
    struct dag_edge_type
    {
        string from_marker;
        string to_marker;
        double weight;
        double confidence; //the confidence in the accuracy of the weight for the edge
    };
    
    struct adj_edge_type { // the satellite data used in graph adjacency list representation
        int to_id;
        double weight;
        double confidence;
    };
    
    class linkage_group_DAG;
    
    /*linakge group in bins representation*/
    class linkage_group_bin
    {
        private:
            double lower_bound;
            double upper_bound;
            int number_of_bins;
            int number_of_markers;
            string description; // used to store the Map name
            double confidence; // the confidence in the accuracy of the LG
            // for bin 0 to the last bin
            // the distance is in CMs
            vector<double> distance_btw_adjacent_pairs; 
            vector<vector<string> > markers_in_bins; 
           
        public:

            linkage_group_bin();
            linkage_group_bin(const linkage_group_bin& _old_lg_bin);
            linkage_group_bin& operator=(const linkage_group_bin& _old_lg_bin);
            
            void initialize(int _in_number_of_bins, 
                            int _in_number_of_markers, 
                            double _in_lower_bound, 
                            double _in_upper_bound,
                            string _description, 
                            double _confidence,
                            const vector<double> & _in_distance_btw_adjacent_pairs, 
                            const vector<vector<string> > & _in_markers_in_bins );
                            
            // false means in the reverse direction and true means in the backward direction. 
            // The return object will be allocated within the function call
            // The callee is responsible for deleting the object after use
            linkage_group_DAG* construct_linkage_group_DAG(bool orientation) const; 
            
            // The following function returns the set of markers that appeared in the LG
            // The callee will allocate the set<string> object inside the function
            // The caller is responsible for deleting the object once it has finished using it
            // Otherwise, memory leak will result
            set<string> * get_mrks_set() const;
            
            
            void remove_markers(vector<string> mrkrs_to_remove);
            
            // The following function will generate a fragment according to the input set
            Fragment* GenFrag(const vector<string>& mrks);
            
            void dump() const;
            ~linkage_group_bin();
            
            void flip_orientation();
            
            int get_number_of_bins() const;
            int get_number_of_markers() const;
            double get_lower_bound() const;
            double get_upper_bound() const;
            const vector<double> & get_distance_btw_adjacent_pairs() const;
            const vector<vector<string> >& get_markers_in_bins() const;
            string get_description() const {
                return description;
            }
            double get_confidence() const {
            	return confidence;
            }
        
        friend bool lgs_intersect(const linkage_group_bin & lg1, 
                                  const linkage_group_bin & lg2);
                                  
        friend pair<double,double> pairwise_tau(const linkage_group_bin & lg1, 
                                                const linkage_group_bin & lg2, 
                                                int & common_markers);

    };

    bool lgs_intersect(const linkage_group_bin & lg1, const linkage_group_bin & lg2);
    pair<double,double> pairwise_tau(const linkage_group_bin & lg1, const linkage_group_bin & lg2, int& common_markers);
    
    class linkage_group_graph_representation;
        
    // A linkage group is a DAG in the adjacency list representation
    // DAG here is a misnomer, the graph may contain cycles if it is obtained by combining multiple dags 
    class linkage_group_DAG 
    {
        private:
            set<string> markers;
            vector<dag_edge_type> edges;
            
        public:
            linkage_group_DAG( );
            void initialize( const set<string>& _in_markers, const vector<dag_edge_type> & _in_edges);
            
            // The returned object is allocated within the function call and needs to be 
            // deallocated by the caller later otherwise memory leak will result
            linkage_group_graph_representation * construct_linkage_group_graph_representation() const;
            
            // The following function will return a reference to the edges member object
            const vector<dag_edge_type> & get_edges() const;
            
            ~linkage_group_DAG();

        friend linkage_group_DAG* combine_lg_DAGs(const linkage_group_DAG* lg_dag1, const linkage_group_DAG* lg_dag2);
    };
    
    // The following function will combine two linkage_group_DAG's to return one single object
    // The two input objects are untouched
    // The callee will allocate the return object inside the function
    // The caller is responsible for deleting the object after use, otherwise memory leak will result
    linkage_group_DAG* combine_lg_DAGs(const linkage_group_DAG* lg_dag1, const linkage_group_DAG* lg_dag2);

    //this graph may contain cycles
    class linkage_group_graph_representation 
    {
        private:
            // The following three data members are input variables
            int number_of_markers;
            vector<string> marker_ids;
            //the second element of a pair is the weight of the edge
            vector<vector<adj_edge_type> > adjacency_list; 
                      
            // The following data members are used by the function that follow
            int num_scc;
            vector<vector<int> > sccs;
            vector<int> singletons;
            // The following functions will be called at the end of initialize function
            void compute_scc();
            
            // the following two data members are to be used by the function that follow
            // true meaning the edge is redundent
            vector<vector<bool> > redundent_edges; 
            int total_number_of_redundent_edges;
            // The following step will remove redundent edges from adjacency_list
            // The function first extract a sub-graph induced by the set of singletons 
            // The sub-graph induced by the set of singletons is guaranteed to be a DAG
            // The function then try to remove redundent edges 
            void minimize_graph(); 
            
            
            // The following data structure and data members are to be used by the function that follow
            struct scc_edge_type{
                int from_id;
            	int to_id;
            	double weight;
            	double confidence;
            	bool redundent;
            };
            vector<vector<scc_edge_type> > scc_adj_list;
            int num_redundent_scc_edge;
            // contract each strongly connected component into one single node
            // generate a DAG in SCC
            void contract_into_scc();
            void minimize_scc();
            
            struct TP_Sort_NS {
                int node_id;
                int in_degree;
            };
            struct compare_NS {
                bool operator()(const TP_Sort_NS& ns1, const TP_Sort_NS& ns2) const {
                    return ns1.in_degree > ns2.in_degree;
                }
            };
            // step 0: initialize 
            void initialize_linearization();
            // step 1: compute pair-wise distances
            void compute_pair_wise_distances();
            // step 2: impute a distance between a pair of markers with no relationship 
            void impute_missing_distance();
            // select a node from a vector a nodes 
            int select_node(const vector<int>& top_nodes);
            // step 3: linearize the dag with modified topological sort
            void linearize_graph();
            vector<int> linear_order; // the linear order by customized topological sort algorithm
            vector<int> tp_linear_order; // the linearized order by ordinary topological sort algorithm
            vector<vector<double> > pair_wise_distances;
                    
        public:
            linkage_group_graph_representation();
            void initialize(int _in_number_of_markers, 
                            const vector<string> & _in_marker_ids, 
                            const vector<vector<adj_edge_type> > _in_adjacency_list);
            int get_num_of_markers() const;
            int get_num_of_scc() const;
            int get_num_of_singletons() const;
            int get_total_number_of_redundent_edges() const;
            pair<string,string> output_dot();
            
            void compute_linear_order();
            
            string output_dot_scc(LG_CLUSTER* lg_clster); // A different method of drawing the graph
            string output_linear_order(); // output the linearized graph in dot file format

            string output_linear_map(); // output the final linearized map in the map chart format
            ~linkage_group_graph_representation();                
    };

    
    class individual_mapping_ppl
    {
        private:
            int number_of_lgs;
            vector<linkage_group_bin> lgs;
            
            // convert recombination fractions into distances in CM

        public:
            individual_mapping_ppl();
            
            /*on error, return -1, otherwise return 0*/
            int read_from_file(string file_path, double map_confidence, string map_name);
            const vector<linkage_group_bin>& get_lgs() const;
            int get_number_of_lgs() const;
            void dump() const;
            ~individual_mapping_ppl();
    };
}// end of namespace

#endif


