/*
 *  condense_markers_into_bins.cpp
 *  consensusmap
 *
 *  Created by yonghui on 5/12/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 
#include "condense_markers_into_bins.h"

using namespace std;
using namespace consensus_map;
using namespace boost;

void consensus_map::consolidate_markers_into_bins(vector<const linkage_group_bin*>& _in_p_to_array_lgs, 
                                                  int _in_number_of_lgs, 
                                                  vector<linkage_group_bin*>& _out_p_to_array_lgs, 
                                                  const linkage_group_DAG & lg_dag)
{

    set<string> universal_mrk_id_str;
    map<string,int> mrk_str_2_id;
    vector<string> mrk_ids;
    int total_number_of_mrks;
    for (int ii = 0; ii<_in_number_of_lgs; ii++ )
    {
        const linkage_group_bin* crt_p_2_lg = _in_p_to_array_lgs[ii];
        const vector< vector<string> > & crt_markers_in_bins = crt_p_2_lg->get_markers_in_bins();
        for (int jj = 0 ; jj < crt_p_2_lg->get_number_of_bins(); jj++)
        {
            const vector<string> & tmp_bin = crt_markers_in_bins[jj];
            universal_mrk_id_str.insert(tmp_bin.begin(), tmp_bin.end());
        }
    }
    int counter = 0;
    for(set<string>::iterator iter1 = universal_mrk_id_str.begin(); iter1 != universal_mrk_id_str.end(); iter1++)
    {
        mrk_str_2_id[*iter1] = counter;
        mrk_ids.push_back(*iter1);
        counter = counter + 1;
    }
    total_number_of_mrks = counter;
    
    
    /*
    step -1: compute the transitive closure
    */
    typedef boost::adjacency_list<vecS,vecS,directedS,no_property,property< edge_weight_t, int> > Graph_distance;
    Graph_distance G_D(total_number_of_mrks);
    const vector<dag_edge_type> & dag_edges = lg_dag.get_edges(); 
    for (vector<dag_edge_type>::const_iterator iter1 = dag_edges.begin(); iter1 != dag_edges.end(); iter1++)
    {
        int from = mrk_str_2_id[iter1->from_marker];
        int to = mrk_str_2_id[iter1->to_marker];
        add_edge(from, to, G_D);
    }
    
    property_map <Graph_distance, edge_weight_t >::type w = get(edge_weight, G_D);
    graph_traits<Graph_distance>::edge_iterator e, e_end;
    for (boost::tie(e, e_end) = boost::edges(G_D); e != e_end; ++e)
    {
        w[*e] = 1;
    }    

    vector<vector<int> > D(total_number_of_mrks, vector<int>(total_number_of_mrks,-1));

    johnson_all_pairs_shortest_paths(G_D, D);    
    
    bool** d = new bool*[total_number_of_mrks];
    for (int ii = 0; ii < total_number_of_mrks; ii++)
    {
        d[ii] = new bool[total_number_of_mrks];
    }
    
    for (int ii =0 ; ii < total_number_of_mrks; ii++)
    {
        for (int jj = 0 ; jj < total_number_of_mrks; jj++)
        {
            if (D[ii][jj] < 0)
            {
                cout << "ERROR! the distance is not as expected" << endl;
            }
            if (D[ii][jj] < total_number_of_mrks +1)
            {
                d[ii][jj] = true;
            }
            else
            {
                d[ii][jj] = false;
            }
        }
    }
    
    /*end of computing the transitive closure for the graph*/
    

    /*
    step 0: break the markers into connected components
    */
    typedef adjacency_list <vecS, vecS, undirectedS> Graph;
    Graph markers_graph(total_number_of_mrks);
        
    for (int ii = 0; ii<_in_number_of_lgs; ii++ )
    {
        const linkage_group_bin* crt_p_2_lg = _in_p_to_array_lgs[ii];
        const vector< vector<string> > & crt_markers_in_bins = crt_p_2_lg->get_markers_in_bins();
        for (int jj = 0 ; jj < crt_p_2_lg->get_number_of_bins(); jj++)
        {
            const vector<string> & tmp_bin = crt_markers_in_bins[jj];
            int size = tmp_bin.size();
            int id1 = mrk_str_2_id[tmp_bin[0]];
            for (int kk = 1; kk < size; kk++)
            {
                int id2 = mrk_str_2_id[tmp_bin[kk]];
                add_edge(id1,id2,markers_graph);
            }
        }
    }
    vector<int> cc_ids(total_number_of_mrks, -1);
    int num = connected_components(markers_graph, &cc_ids[0]);
    vector<vector<int > > connected_components(num);
    
    for (int ii = 0; ii < total_number_of_mrks; ii++ )
    {
        connected_components[cc_ids[ii]].push_back(ii);
    }
    
    /*
    step 1: Condense the markers into bins.
    */
    
    vector< vector<int> > co_segregating_bins;
    for (int ii = 0; ii < num; ii++)
    {
        /*
        For each connected component, break it into co-segregating bins. 
        */
        set<int> remaining_nodes;
        remaining_nodes.insert(connected_components[ii].begin(), connected_components[ii].end());
        vector<int> crt_bin;
        while( remaining_nodes.size() != 0 )
        {
            crt_bin.clear();
            for (set<int>::iterator iter1 = remaining_nodes.begin(); iter1 != remaining_nodes.end(); iter1++)
            {
                bool compatible = true;
                for (vector<int>::iterator iter2 = crt_bin.begin(); iter2 != crt_bin.end(); iter2++)
                {
                    /* check if *iter1 is compatible with *iter2 */
                    if((d[*iter1][*iter2]) or (d[*iter2][*iter1]))
                    {
                        compatible = false;
                        break;
                    }
                }
                if (compatible)
                {
                    /*Then, we need to modify the dag to reflect the fact that *iter1 has been added to the super ndoe*/
                    set<pair<int,int> > edges_to_be_added;
                    
                    if( crt_bin.size() != 0 )
                    {
                        set<int> old_parents;
                        set<int> new_parents;
                        set<int> old_children;
                        set<int> new_children;
                        int old_mrk = crt_bin[0];
                        int new_mrk = *iter1;
                        for (int jj5 = 0; jj5 < total_number_of_mrks; jj5++)
                        {
                            if (d[jj5][old_mrk]) old_parents.insert(jj5);
                            if (d[jj5][new_mrk]) new_parents.insert(jj5);
                            if (d[old_mrk][jj5]) old_children.insert(jj5);
                            if (d[new_mrk][jj5]) new_children.insert(jj5);
                        }
                        set<int> additional_new_parents = new_parents;
                        set<int> additional_old_parents = old_parents;
                        for (set<int>::iterator iter5 = old_parents.begin(); iter5 != old_parents.end(); iter5++)
                        {
                            additional_new_parents.erase(*iter5);
                        }
                        for (set<int>::iterator iter5 = new_parents.begin(); iter5 != new_parents.end(); iter5++)
                        {
                            additional_old_parents.erase(*iter5);
                        }
                        for (set<int>::iterator iter6 = additional_new_parents.begin(); iter6 != additional_new_parents.end(); iter6++)
                        {
                            for (set<int>::iterator iter7 = old_children.begin(); iter7 != old_children.end(); iter7++)
                            {
                                edges_to_be_added.insert(make_pair(*iter6, *iter7));
                            }
                            for (vector<int>::iterator iter7 = crt_bin.begin(); iter7 != crt_bin.end()  ; iter7++)
                            {
                                edges_to_be_added.insert(make_pair(*iter6, *iter7));
                            }
                        }
                        for (set<int>::iterator iter6 = additional_old_parents.begin(); iter6 != additional_old_parents.end(); iter6++)
                        {
                            for (set<int>::iterator iter7 = new_children.begin(); iter7 != new_children.end(); iter7++)
                            {
                                edges_to_be_added.insert(make_pair(*iter6, *iter7));
                            }
                            edges_to_be_added.insert(make_pair(*iter6, new_mrk));
                        }
                        for (set<int>::iterator iter7 = old_children.begin(); iter7 != old_children.end(); iter7++)
                        {
                            edges_to_be_added.insert(make_pair(new_mrk, *iter7));
                        }
                        for (vector<int>::iterator iter7 = crt_bin.begin(); iter7 != crt_bin.end(); iter7++)
                        {
                            for (set<int>::iterator iter8 = new_children.begin(); iter8 != new_children.end(); iter8++ )
                            {
                                edges_to_be_added.insert(make_pair(*iter7, *iter8));
                            }
                        }
                    }
                    for (set<pair<int,int> >::iterator iter8 = edges_to_be_added.begin(); iter8 != edges_to_be_added.end(); iter8++)
                    {
                        d[iter8->first][iter8->second] = true;
                    }
                    
                    crt_bin.push_back(*iter1);
                    
                }
            }
            co_segregating_bins.push_back(crt_bin);
            for (vector<int>::iterator iter1 = crt_bin.begin(); iter1 != crt_bin.end(); iter1++)
            {
                remaining_nodes.erase(*iter1);
            }
            
        }

    }
    
    /*
    step 2: construct a map relates original id to the super node id
    */
    map<int,int> id_2_super_id;
    for (int ii1 = 0 ; ii1 < co_segregating_bins.size(); ii1++)
    {
        for (vector<int>::iterator iter1 = co_segregating_bins[ii1].begin();  iter1 != co_segregating_bins[ii1].end(); iter1++)
        {
            id_2_super_id[*iter1] = ii1;
        }
    }
    
    
    vector<string> super_node_str_ids;
    for (int ii1 = 0 ; ii1 < co_segregating_bins.size(); ii1++)
    {
        string tmp = mrk_ids[co_segregating_bins[ii1][0]];
        for (int ii2 = 1; ii2 < co_segregating_bins[ii1].size(); ii2++)
        {
            tmp = tmp + ',' + mrk_ids[co_segregating_bins[ii1][ii2]];
        }
        super_node_str_ids.push_back(tmp);
    }
    
    /*
    reconstruct the lgs in the super-node representation
    */
    for (int ii3 = 0 ; ii3 < _in_number_of_lgs; ii3++)
    {
        const vector<vector<string> > & original_bins = _in_p_to_array_lgs[ii3]->get_markers_in_bins();
        vector<vector<string> > new_bins;
        int number_of_markers = 0;
        for (int ii4 = 0; ii4 < original_bins.size(); ii4++)
        {
            set<int> super_node_ids;
            for (int ii5 = 0 ; ii5 < original_bins[ii4].size(); ii5++)
            {
                int tmp_id = id_2_super_id[mrk_str_2_id[original_bins[ii4][ii5]]];
                super_node_ids.insert(tmp_id);
            }
            vector<string> crt_bin;
            for (set<int>::iterator iter1 = super_node_ids.begin(); iter1 != super_node_ids.end(); iter1++)
            {
                crt_bin.push_back(super_node_str_ids[*iter1]);
            }
            new_bins.push_back(crt_bin);
            number_of_markers = number_of_markers + crt_bin.size();
        }

        int number_of_bins = _in_p_to_array_lgs[ii3]->get_number_of_bins();
        double lower_bound = _in_p_to_array_lgs[ii3]->get_lower_bound();
        double upper_bound = _in_p_to_array_lgs[ii3]->get_upper_bound();
        string description = _in_p_to_array_lgs[ii3]->get_description();
        double confidence = _in_p_to_array_lgs[ii3]->get_confidence();
        const vector<double> & distance_btw_adjacent_pairs = _in_p_to_array_lgs[ii3]->get_distance_btw_adjacent_pairs();
        linkage_group_bin * new_lg = new linkage_group_bin();
        new_lg->initialize(number_of_bins, 
                           number_of_markers, 
                           lower_bound, 
                           upper_bound, 
                           description,
                           confidence,
                           distance_btw_adjacent_pairs, 
                           new_bins);
        _out_p_to_array_lgs[ii3] = new_lg;
    }
    
    /*
    deallocate memory
    */    
    for (int ii = 0; ii < total_number_of_mrks; ii++)
    {
        delete[] d[ii];
    }
    delete[] d;
    
    return;
};
