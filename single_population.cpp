/*
 *  single_population.cpp
 *  consensusmap
 *
 *  Created by yonghui on 5/12/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include <cassert>
#include "single_population.h"

namespace consensus_map{

    linkage_group_bin::linkage_group_bin()
    {
        lower_bound = 0.0;
        upper_bound = 0.0;
        number_of_bins = 0;
        number_of_markers = 0;
        description = "";
        confidence = 0.0;
        distance_btw_adjacent_pairs.clear(); //for bin 0 to the last bin
        markers_in_bins.clear(); 
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    linkage_group_bin::linkage_group_bin(const linkage_group_bin& _old_lg_bin){
        (*this) = _old_lg_bin;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    linkage_group_bin& linkage_group_bin::operator=(const linkage_group_bin& _old_lg_bin){
        lower_bound = _old_lg_bin.get_lower_bound();
        upper_bound = _old_lg_bin.get_upper_bound();
        number_of_bins = _old_lg_bin.get_number_of_bins();
        number_of_markers = _old_lg_bin.get_number_of_markers();
        distance_btw_adjacent_pairs = _old_lg_bin.get_distance_btw_adjacent_pairs(); 
        markers_in_bins = _old_lg_bin.get_markers_in_bins(); 
        description = _old_lg_bin.get_description();
        confidence = _old_lg_bin.get_confidence();
        return (*this);
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    void linkage_group_bin::remove_markers(vector<string> mrkrs_to_remove){
        
        set<string> set_mrkrs_to_remove;
        set_mrkrs_to_remove.insert(mrkrs_to_remove.begin(), mrkrs_to_remove.end());
        
        // assert mrkrs_to_remove is a singleton set
        assert(set_mrkrs_to_remove.size() == mrkrs_to_remove.size());
        
        int remove_count = 0;
        
        map<string,int> mrk_2_level;
        set<int> old_levels;
        map<int, int> old_2_new_level;
        
        assert(number_of_bins == markers_in_bins.size());
        assert(number_of_bins == distance_btw_adjacent_pairs.size() + 1);
        for (int ii = 0; ii < number_of_bins; ii++) {
            for (vector<string>::iterator iter1 = markers_in_bins[ii].begin(); 
                iter1 != markers_in_bins[ii].end();
                ++iter1) {
                if (set_mrkrs_to_remove.find(*iter1) != set_mrkrs_to_remove.end()) {
                    remove_count++;
                } else {
                    mrk_2_level[*iter1] = ii;
                    old_levels.insert(ii);                 
                }
            }
        }
         
        // assert the number of markers to be removed is equal to the input size
        assert(remove_count == mrkrs_to_remove.size());
        
        vector<int> vec_old_levels;
        vec_old_levels.insert(vec_old_levels.begin(), old_levels.begin(), old_levels.end());
        sort(vec_old_levels.begin(), vec_old_levels.end());
        for (int ii = 0; ii < vec_old_levels.size(); ii++) {
            old_2_new_level[vec_old_levels[ii]] = ii;
        }
                
        vector<vector<string> > new_markers_bins(old_levels.size());
        for (map<string, int>::const_iterator iter2 = mrk_2_level.begin(); iter2 != mrk_2_level.end(); ++iter2) {
                string marker = iter2->first;
                int new_level_id = old_2_new_level[iter2->second];
                assert(new_level_id < old_levels.size());
                new_markers_bins[new_level_id].push_back(marker);
        }        
        
        vector<double> new_adj_distances;
        if (old_levels.size() > 1) {
            new_adj_distances.resize(old_levels.size() -1);
            for (int ii = 0; ii < new_adj_distances.size(); ii++) {
                new_adj_distances[ii] = 0;
                for (int jj = vec_old_levels[ii]; jj < vec_old_levels[ii + 1]; jj++) {
                    new_adj_distances[ii] = new_adj_distances[ii] + distance_btw_adjacent_pairs[jj];
                }
            }
        }
        number_of_bins = old_levels.size();
        number_of_markers = number_of_markers - mrkrs_to_remove.size();
        markers_in_bins = new_markers_bins;
        distance_btw_adjacent_pairs = new_adj_distances;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
    void linkage_group_bin::flip_orientation(){
        int num_dists = distance_btw_adjacent_pairs.size();
        vector<double> tmp_dist = distance_btw_adjacent_pairs;
        int num_bins = markers_in_bins.size();
        vector<vector<string> > tmp_markers = markers_in_bins;
        for (int i = num_bins - 1; i >= 0; i--) {
            markers_in_bins[num_bins - i - 1] = tmp_markers[i];
        }
        for (int i = num_dists - 1; i >= 0; i--) {
            distance_btw_adjacent_pairs[num_dists - i - 1] = tmp_dist[i];
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    linkage_group_bin::~linkage_group_bin()
    {
        return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    linkage_group_DAG* linkage_group_bin::construct_linkage_group_DAG(bool orientation) const
    {
            vector<dag_edge_type> edges;
            set<string> markers;
            
            /*to insert the markers in the first bin into "markers" object*/
            for (vector<string>::const_iterator iter1 = (markers_in_bins[0]).begin(); 
                 iter1 != (markers_in_bins[0]).end(); 
                 iter1++)
            {
                markers.insert(*iter1);
            }
            
            dag_edge_type one_edge;
            
            /*
            The following nested for loop will insert no markers into the "markers" object if there is only one bin
            */
            for (int ii = 0; ii < number_of_bins-1; ii++)
            {
                for (vector<string>::const_iterator iter1 = (markers_in_bins[ii]).begin(); 
                     iter1 != (markers_in_bins[ii]).end(); 
                     iter1++)
                {
                    markers.insert(*iter1);
                    for ( vector<string>::const_iterator iter2 = (markers_in_bins[ii+1]).begin(); 
                          iter2 != (markers_in_bins[ii+1]).end(); 
                          iter2++)
                    {
                        markers.insert(*iter2);
                        one_edge.weight = distance_btw_adjacent_pairs[ii];  
                        one_edge.confidence = confidence;
                        if (orientation) // in the forward direction
                        {
                            one_edge.from_marker = *iter1;
                            one_edge.to_marker = *iter2;
                        }
                        else 
                        {
                            one_edge.from_marker = *iter2;
                            one_edge.to_marker = *iter1;
                        }
                        edges.push_back(one_edge);
                    }
                }                
            }

            linkage_group_DAG * to_return = new linkage_group_DAG;
            to_return->initialize(markers, edges);
            return to_return;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    void linkage_group_bin::initialize(int _in_number_of_bins, 
                                       int _in_number_of_markers, 
                                       double _in_lower_bound, 
                                       double _in_upper_bound, 
                                       string _description,
                                       double _confidence,
                                       const vector<double> & _in_distance_btw_adjacent_pairs, 
                                       const vector<vector<string> > & _in_markers_in_bins )
    {
        lower_bound = _in_lower_bound;
        upper_bound = _in_upper_bound;
        number_of_bins = _in_number_of_bins;
        number_of_markers = _in_number_of_markers;
        description = _description;
        confidence = _confidence;
        distance_btw_adjacent_pairs = _in_distance_btw_adjacent_pairs;
        markers_in_bins = _in_markers_in_bins;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    set<string> * linkage_group_bin::get_mrks_set() const
    {
        set<string>* to_return = new set<string>;
        for (int ii = 0 ; ii < number_of_bins; ii++)
        {
            to_return->insert(markers_in_bins[ii].begin(),markers_in_bins[ii].end());
        }
        return to_return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    Fragment* linkage_group_bin::GenFrag(const vector<string>& mrks){
        map<string,int> mrk_2_level;
        assert(number_of_bins == markers_in_bins.size());
        assert(number_of_bins == distance_btw_adjacent_pairs.size() + 1);
        for (int ii = 0; ii < number_of_bins; ii++) {
            for (vector<string>::iterator iter1 = markers_in_bins[ii].begin(); 
                 iter1 != markers_in_bins[ii].end();
                 ++iter1) {
                 mrk_2_level[*iter1] = ii;
            }
        }
        
        set<int> related_levels;
        for (vector<string>::const_iterator iter2 = mrks.begin(); iter2 != mrks.end(); ++iter2) {
            if ( mrk_2_level.find(*iter2) != mrk_2_level.end()) {
                related_levels.insert(mrk_2_level[*iter2]);
            }
        }
        
        if (related_levels.size() == 0) return NULL;
        
        int max_level = 0;
        int min_level = number_of_bins;
        
        for (set<int>::iterator iter3 = related_levels.begin();
             iter3 != related_levels.end();
             ++iter3) {
             if (max_level < *iter3) max_level = *iter3;
             if (min_level > *iter3) min_level = *iter3;
        }
        
        // it is expected that a fragment constitutes of mrks from consecutive levels 
        assert((max_level - min_level) < related_levels.size());
        
        int num_int_mrks = 0;
        vector<vector<string> > frag_bins(related_levels.size()); 
        for (vector<string>::const_iterator iter2 = mrks.begin(); iter2 != mrks.end(); ++iter2) {
            if ( mrk_2_level.find(*iter2) != mrk_2_level.end()) {
                frag_bins[mrk_2_level[*iter2] - min_level].push_back(*iter2);
                if ((mrk_2_level[*iter2] != min_level) and (mrk_2_level[*iter2] != max_level)) {
                    num_int_mrks = num_int_mrks + 1;
                }
            }
        }        
        
        int num_int_mrks_2 = 0;
        for (int ii = min_level + 1; ii <= max_level - 1; ii++) {
            num_int_mrks_2 = num_int_mrks_2 + markers_in_bins[ii].size();
        }
        
        // it is expected that a fragment will contain all markers from intermediate levels
        assert(num_int_mrks == num_int_mrks_2);
        
        vector<double> frag_dists;
        for (int ii = min_level; ii < max_level; ii++) {
            frag_dists.push_back(distance_btw_adjacent_pairs[ii]);
        }
        
        Fragment* to_return = new Fragment(frag_bins, frag_dists, confidence, description);
        return to_return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int linkage_group_bin::get_number_of_bins() const
    {
        return number_of_bins;
    };
    int linkage_group_bin::get_number_of_markers() const
    {
        return number_of_markers;
    };
    double linkage_group_bin::get_lower_bound() const
    {
        return lower_bound;
    };
    double linkage_group_bin::get_upper_bound() const
    {
        return upper_bound;
    };
    const vector<double> & linkage_group_bin::get_distance_btw_adjacent_pairs() const
    {
        return distance_btw_adjacent_pairs;
    };
    const vector<vector<string> >& linkage_group_bin::get_markers_in_bins() const
    {
        return markers_in_bins;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void linkage_group_bin::dump() const
    {
        cout << "==============================================" << endl;
        cout << "upperbound:" << upper_bound << ", lowerbound:" << lower_bound  
             << ", confidence:" << confidence << endl
             << "description: " << description << endl
             << ", number_of_bins:" << number_of_bins 
             << ", number_of_markers:" << number_of_markers << endl;
             
        cout << "the pairwise distances:" << endl;
        for (int ii = 0 ; ii < number_of_bins -1; ii++)
        {
            cout << distance_btw_adjacent_pairs[ii] << '\t' ;
        }
        cout << endl;
        for (int ii = 0 ; ii < number_of_bins; ii++)
        {
            for (vector<string>::const_iterator iter1 = markers_in_bins[ii].begin(); 
                 iter1 != markers_in_bins[ii].end(); 
                 iter1++)
            {
                cout << *iter1 << '\t' ;
            }
            cout << endl;
        }
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    bool lgs_intersect(const linkage_group_bin & lg1, const linkage_group_bin & lg2)
    {
        set<string> mrks_set1;
        set<string> mrks_set2;
        for (int ii =0 ; ii < lg1.number_of_bins; ii++)
        {
            mrks_set1.insert((lg1.markers_in_bins[ii]).begin(),(lg1.markers_in_bins[ii]).end() );
        }
        for (int ii =0; ii < lg2.number_of_bins; ii++)
        {
            mrks_set2.insert((lg2.markers_in_bins[ii]).begin(),(lg2.markers_in_bins[ii]).end() );
        }
        bool found = false;
        set<string>::iterator end_of_set2 = mrks_set2.end();
        int num_common = 0;
        set<string> common_markers;
        for (set<string>::iterator iter1=mrks_set1.begin(); iter1!= mrks_set1.end(); iter1++)
        {
            if ( mrks_set2.find(*iter1) != end_of_set2)
            {
                found = true;
		num_common = num_common + 1;
		common_markers.insert(*iter1);
            }
        }
        // edited by yonghui on Feb 20
        if ((num_common > 0) and (num_common <= klgMinCommon)) {
            cout << "caution: two lgs share too few markers in common:" ;
            for (set<string>::iterator iter1=common_markers.begin(); iter1!= common_markers.end(); iter1++) {
                cout << *iter1 << '\t';
            }
            cout << endl;
        }
        return found;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    // The following function will return tau between two LGs
    // The first element of the returned pair is the tau if the two LGs are compared as they are
    // The second element of the returned pair is the tau when the direction of one of the LGs is flipped 
    pair<double,double> pairwise_tau(const linkage_group_bin & lg1, 
                                     const linkage_group_bin & lg2, 
                                     int& common_markers)
    {
        vector<int> rank1;
        vector<int> rank2;
        
        /*1. get the set of overlapping markers*/
        map<string,int>  map1;

        for (int ii = 0 ; ii < lg1.number_of_bins; ii++)
        {
            int number_of_markers_in_bin = ((lg1.markers_in_bins)[ii]).size();
            for (int jj = 0 ; jj < number_of_markers_in_bin ; jj++)
            {
                map1[lg1.markers_in_bins[ii][jj]] = ii;
            }
        }

        map<string,int>::iterator map1_end = map1.end();
        for (int ii = 0 ; ii < lg2.number_of_bins; ii++)
        {
            int number_of_markers_in_bin = lg2.markers_in_bins[ii].size();
            for (int jj = 0; jj < number_of_markers_in_bin ; jj++)
            {
                map<string,int>::iterator tmp_iter = map1.find(lg2.markers_in_bins[ii][jj]);
                if (not (tmp_iter == map1_end))
                {
                    rank1.push_back(tmp_iter->second);
                    rank2.push_back(ii);
                }
            }
        }
        common_markers = rank1.size();
        /*2. calculate the pair-wise tau*/
        int number_of_common_markers = rank1.size();
        if (number_of_common_markers <= 1)
        {
            return make_pair(0,0);
        }
        else
        {
            double concordant_pairs = 0 ;
            double conflict_pairs = 0 ;
            for (int ii = 0 ; ii < number_of_common_markers; ii++)
            {
                for (int jj = ii+1; jj < number_of_common_markers; jj++)
                {
                    if ((rank1[ii]-rank1[jj]) * (rank2[ii] - rank2[jj]) >=0)
                    {
                        concordant_pairs = concordant_pairs + 1;
                    }
                    if ((rank1[ii]-rank1[jj]) * (rank2[ii] - rank2[jj]) <=0)
                    {
                        conflict_pairs = conflict_pairs + 1;
                    }
                }
            }
            double forward_tau = 4*concordant_pairs/(number_of_common_markers * (number_of_common_markers -1.0)) - 1.0;
            double backward_tau = 4*conflict_pairs/(number_of_common_markers * (number_of_common_markers -1.0)) - 1.0;
            return make_pair(forward_tau,backward_tau);
        }

    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    linkage_group_DAG* combine_lg_DAGs(const linkage_group_DAG* lg_dag1, const linkage_group_DAG* lg_dag2) {
        const set<string> & mark_set1 = lg_dag1->markers;
        const set<string> & mark_set2 = lg_dag2->markers;
        const vector<dag_edge_type>& edges1 = lg_dag1->edges;
        const vector<dag_edge_type>& edges2 = lg_dag2->edges;
        
        set<string> markers_combined(mark_set1);
        markers_combined.insert(mark_set2.begin(), mark_set2.end());
        
        /*Now combine the edges set*/
        map< pair<string,string>, pair<double,double> > tmp_map;
        for (vector<dag_edge_type>::const_iterator iter1 = edges1.begin(); iter1 != edges1.end(); iter1++ )
        {
            tmp_map[make_pair(iter1->from_marker,iter1->to_marker)] = make_pair(iter1->weight,iter1->confidence);
        }
        
        for (vector<dag_edge_type>::const_iterator iter2 = edges2.begin(); iter2 != edges2.end(); iter2++ )
        {
            string from_marker= iter2->from_marker;
            string to_marker = iter2->to_marker;
            double weight = iter2->weight;
            double confidence = iter2 -> confidence; 
            if (tmp_map.find(make_pair(from_marker,to_marker)) == tmp_map.end())
            {
                tmp_map[make_pair(from_marker, to_marker)] = make_pair(weight, confidence);
            }
            else 
            {
                double old_confidence = tmp_map[make_pair(from_marker, to_marker)].second;
                double old_weight = tmp_map[make_pair(from_marker, to_marker)].first;
                double new_confidence = old_confidence + confidence;
                double new_weight = ( old_weight * old_confidence + weight * confidence ) / new_confidence ; 
                tmp_map[make_pair(from_marker, to_marker)] = make_pair(new_weight, new_confidence);
            }
        }
        vector<dag_edge_type> edges_combined;
        for (map<pair<string,string>, pair<double, double> >::iterator iter1 = tmp_map.begin(); 
             iter1 != tmp_map.end(); 
             iter1++)
        {
            string from_marker = (iter1->first).first ;
            string to_marker = (iter1->first).second ;
            double weight = (iter1->second).first ;
            double confidence = (iter1->second).second ; 
            dag_edge_type tmp_edge = {from_marker, to_marker, weight, confidence};
            edges_combined.push_back(tmp_edge);
        }
        linkage_group_DAG* to_return = new linkage_group_DAG();
        to_return->initialize(markers_combined, edges_combined);
        return to_return;
        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    linkage_group_DAG::linkage_group_DAG( )
    {
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    void linkage_group_DAG::initialize( const set<string>& _in_markers, const vector<dag_edge_type> & _in_edges) 
    {
        markers = _in_markers;
        edges = _in_edges;
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    linkage_group_graph_representation * linkage_group_DAG::construct_linkage_group_graph_representation() const
    {
        vector<string> markers_set;
        map<string,int> marker_id;
        int counter = 0 ;
        for (set<string>::const_iterator iter1 = markers.begin(); iter1 != markers.end(); iter1++ )
        {
            markers_set.push_back(*iter1);
            marker_id[*iter1] = counter;
            counter = counter + 1;
        }
        vector<vector<adj_edge_type> > adjacency_list(counter);
        for (vector<dag_edge_type>::const_iterator iter1 = edges.begin(); iter1 != edges.end(); iter1++)
        {
            int from_id =  marker_id[iter1->from_marker];
            int to_id =  marker_id[iter1->to_marker];
            double weight = iter1->weight;
            double confidence = iter1->confidence;
            adj_edge_type tmp_edge = {to_id,weight,confidence};
            adjacency_list[from_id].push_back(tmp_edge);
        }
        linkage_group_graph_representation * to_return = new linkage_group_graph_representation;
        to_return->initialize(counter, markers_set, adjacency_list);
        return to_return;
    };


    const vector<dag_edge_type> & linkage_group_DAG::get_edges() const
    {
        return edges;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    linkage_group_DAG::~linkage_group_DAG()
    {
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    individual_mapping_ppl::individual_mapping_ppl()
    {
        number_of_lgs = 0;
        lgs.clear();
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    int individual_mapping_ppl::read_from_file(string file_path, double map_confidence, string map_name)
    {
        /*open the file*/
        ifstream inputfile(file_path.c_str());

        vector<int> mrks_in_lgs;
        vector<int> bins_in_lgs;
        vector<double> upperbounds;
        vector<double> lowerbounds;
        vector<vector<double> > distances_btw_adjacent_pairs;
        vector<vector<vector<string> > > lg_bin_mrks;

        /*read the entire file into a string object and then close the file*/
        inputfile.seekg(0, ios_base::end);
        int file_size = inputfile.tellg();
        inputfile.seekg(0, ios_base::beg);
        char * buffer = new char[file_size];
        inputfile.read(buffer, file_size);
        string file_content(buffer);
        delete[] buffer;
        inputfile.close();
        
        number_of_lgs = 0;
        int lg_s = file_content.find(start_lg);
        int lg_e = file_content.find(end_lg);
        while ((lg_s != string::npos) and (lg_e != string::npos)) {
            lg_bin_mrks.push_back(vector<vector<string> >()); // add a new linkage group 
            distances_btw_adjacent_pairs.push_back(vector<double>());
            
            number_of_lgs = number_of_lgs + 1;
            
            string lg_substr = file_content.substr(lg_s + start_lg.length(), lg_e - lg_s - start_lg.length());
            istringstream mrk_bin_ss(lg_substr, ios_base::in);
            
            string marker_id;
            string dist_str;

            string previous_dist = "";
            int total_mrks_ii = 0;
            int total_bins_ii = 0;
            mrk_bin_ss >> marker_id;
            mrk_bin_ss >> dist_str;
            
            while (!mrk_bin_ss.fail())
            {
                total_mrks_ii = total_mrks_ii + 1;
                if (dist_str != previous_dist) { 
                    total_bins_ii = total_bins_ii + 1;
                    lg_bin_mrks[number_of_lgs - 1].push_back(vector<string>()); // create a new bin
                    if (lg_bin_mrks[number_of_lgs - 1].size() != 1) {// not the first bin
                        double dist1 = atof(previous_dist.c_str());
                        double dist2 = atof(dist_str.c_str());
                        assert(dist2 > dist1);
                        distances_btw_adjacent_pairs[number_of_lgs - 1].push_back(dist2 - dist1);
                    }
                    previous_dist = dist_str;
                }
                int number_of_bins = lg_bin_mrks[number_of_lgs - 1].size();
                (lg_bin_mrks[number_of_lgs - 1][number_of_bins - 1]).push_back(marker_id);
                mrk_bin_ss >> marker_id;
                mrk_bin_ss >> dist_str;
            }
            assert(lg_bin_mrks[number_of_lgs - 1].size() == distances_btw_adjacent_pairs[number_of_lgs - 1].size() + 1);
            mrks_in_lgs.push_back(total_mrks_ii);
            bins_in_lgs.push_back(total_bins_ii);
            
            lg_s = file_content.find(start_lg, lg_e + 1);
            lg_e = file_content.find(end_lg, lg_e + 1);            
        }

        upperbounds.resize(number_of_lgs);
        lowerbounds.resize(number_of_lgs);
        for (int ii = 0; ii < number_of_lgs; ii++) { // put dummy upper and lower bound
            lowerbounds[ii] = 0.0;
            upperbounds[ii] = 0.0;
        }

        lgs.resize(number_of_lgs);
        for (int ii = 0 ; ii < number_of_lgs ; ii++)
        {
            char lg_descript[100];
            sprintf(lg_descript, "%s_%d", map_name.c_str(), ii);
            lgs[ii].initialize(bins_in_lgs[ii], 
                               mrks_in_lgs[ii],
                               lowerbounds[ii],
                               upperbounds[ii],
                               lg_descript,
                               map_confidence,
                               distances_btw_adjacent_pairs[ii],
                               lg_bin_mrks[ii]);
        }

        /*if successful, return 0*/
        return 0 ;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    const vector<linkage_group_bin>& individual_mapping_ppl::get_lgs() const
    {
        return lgs;
    };


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    int individual_mapping_ppl::get_number_of_lgs() const
    {
        return number_of_lgs;
    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void individual_mapping_ppl::dump() const
    {
        for (int ii = 0 ; ii < number_of_lgs; ii++)
        {
            lgs[ii].dump();
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    individual_mapping_ppl::~individual_mapping_ppl()
    {
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
    linkage_group_graph_representation::linkage_group_graph_representation()
    {
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void linkage_group_graph_representation::initialize(int _in_number_of_markers, 
                                                        const vector<string> & _in_marker_ids, 
                                                        const vector<vector<adj_edge_type> > _in_adjacency_list)
    {
        number_of_markers = _in_number_of_markers;
        marker_ids = _in_marker_ids;
        adjacency_list = _in_adjacency_list;
        compute_scc();
        minimize_graph();
        contract_into_scc();
        minimize_scc();
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void linkage_group_graph_representation::compute_scc()
    {
        /*clear up*/
        num_scc = 0 ;
        sccs.clear();
        singletons.clear();
            
        typedef boost::adjacency_list<vecS,vecS,directedS> Graph;
        Graph G(number_of_markers);
        for (int ii = 0 ; ii < number_of_markers; ii++) {
            for (vector<adj_edge_type>::iterator iter1 = adjacency_list[ii].begin(); 
                 iter1 != adjacency_list[ii].end(); 
                 iter1++) {
                add_edge(ii,iter1->to_id,G);
            }
        }
        vector<int> components(number_of_markers,-1);
        int num = strong_components(G,&components[0]);
        num_scc = num;
        sccs.resize(num);
        for (int ii = 0; ii < number_of_markers; ii++)
        {
            if (components[ii] >= num)
            {
                cout << "ERROR! strong component id is unexpected" << endl;
                
            }
            else
            {
                sccs[components[ii]].push_back(ii);
            }
        }
        singletons.clear();
        for (int ii = 0; ii < num; ii++)
        {
            if (sccs[ii].size() == 1)
            {
                singletons.push_back(sccs[ii][0]);
            }
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    int linkage_group_graph_representation::get_num_of_markers() const
    {
        return number_of_markers;
    };
    int linkage_group_graph_representation::get_num_of_scc() const
    {
        return num_scc;
    };
    int linkage_group_graph_representation::get_num_of_singletons() const
    {
        return singletons.size();
    };

    int linkage_group_graph_representation::get_total_number_of_redundent_edges() const
    {
        return total_number_of_redundent_edges;
    };

    void linkage_group_graph_representation::minimize_graph()
    {
        /*In the following, the edge weights are assumed to 1, namely the graph is unweighted*/
        typedef boost::adjacency_list<vecS,vecS,directedS,no_property,property< edge_weight_t, int> > Graph;
        Graph sub_G(singletons.size());
        map<int,int> ori_2_new;
        for (int ii = 0; ii < singletons.size(); ii++)
        {
            ori_2_new[singletons[ii]] = ii;
        }
        map<int,int>::iterator end_of_map=ori_2_new.end();
        for (int ii = 0; ii < singletons.size(); ii++)
        {
            for (vector<adj_edge_type>::iterator iter1 = adjacency_list[singletons[ii]].begin(); 
                 iter1 != adjacency_list[singletons[ii]].end(); 
                 iter1 ++)
            {
                if (ori_2_new.find(iter1->to_id) != end_of_map)
                {
                    add_edge(ii,ori_2_new[iter1->to_id],sub_G);
                }
            }
        }
        
        property_map <Graph, edge_weight_t >::type w = get(edge_weight, sub_G);
        graph_traits<Graph>::edge_iterator e, e_end;
        for (boost::tie(e, e_end) = edges(sub_G); e != e_end; ++e)
        {
            w[*e] = 1;
        }
        
        /*Run johnson's all pairs shortest-path algorithm to figure out the transitive closure*/
        int V_size = singletons.size();
        vector<vector<int> > D(V_size, vector<int>(V_size,-1));

        johnson_all_pairs_shortest_paths(sub_G, D);
        
        bool** d = new bool*[V_size];
        for (int ii = 0; ii < V_size; ii++)
        {
            d[ii] = new bool[V_size];
        }

        for (int ii =0 ; ii < V_size; ii++)
        {
            for (int jj = 0 ; jj < V_size; jj++)
            {
                if (D[ii][jj] < 0)
                {
                    cout << "ERROR! the distance is not as expected" << endl;
                }
                if (D[ii][jj] < V_size +1)
                {
                    d[ii][jj] = true;
                }
                else
                {
                    d[ii][jj] = false;
                }
            }
        }
        
        for (int ii = 0 ; ii < V_size; ii++)
        {
            d[ii][ii] = false;
        }
        
        vector<vector<int> > two_path(V_size,vector<int>(V_size,0));
        for (int ii = 0 ; ii < V_size; ii++)
        {
            for (int jj = 0; jj < V_size; jj++)
            {
                for (int kk = 0; kk < V_size; kk++)
                {
                    if (d[ii][kk] and (d[kk][jj]))     
                    {
                        two_path[ii][jj] = 1;
                        break;
                    }
                }
            }
        }
        
        /*delete the redundent edges*/
        total_number_of_redundent_edges = 0 ; 
        redundent_edges.clear();
        redundent_edges.resize(number_of_markers);
        for (int ii = 0; ii < number_of_markers; ii ++) {
            for (vector<adj_edge_type>::iterator iter1 = adjacency_list[ii].begin(); 
                 iter1 != adjacency_list[ii].end(); 
                 iter1++) {
                if ((ori_2_new.find(ii) == end_of_map) or (ori_2_new.find(iter1->to_id) == end_of_map)) {
                    redundent_edges[ii].push_back(false);
                }
                else {
                    int new_ii = ori_2_new[ii];
                    int new_iter1 = ori_2_new[iter1->to_id];
                    if (two_path[new_ii][new_iter1] == 1) {
                        redundent_edges[ii].push_back(true);
                        total_number_of_redundent_edges = total_number_of_redundent_edges + 1;
                    }
                    else {
                        redundent_edges[ii].push_back(false);    
                    }
                }
            }
        }
        
        /*clean up*/
        for (int ii = 0 ; ii < V_size ;ii++)
        {
            delete[] d[ii];
        }
        delete[] d;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    pair<string,string> linkage_group_graph_representation::output_dot()
    {

        ostringstream ss_out1(ios_base::out);
        ostringstream ss_out2(ios_base::out);
        ss_out1 << "digraph G {" << endl;
        ss_out2 << "digraph G {" << endl;
        for (int ii = 0; ii < number_of_markers; ii++)
        {
            for (int jj = 0 ; jj < adjacency_list[ii].size(); jj++)
            {
                if (not redundent_edges[ii][jj])
                {
                    double weight = adjacency_list[ii][jj].weight; 
                    char buffer[100];
                    sprintf(buffer, "%.2f",weight);
                    ss_out1 << "\"" << marker_ids[ii] << "\"" 
                            << " -> " 
                            << "\"" << marker_ids[adjacency_list[ii][jj].to_id] << "\"" 
                            << "[label =" << buffer << "]" << ';' << endl;
                            
                    ss_out2 << "\"" << marker_ids[adjacency_list[ii][jj].to_id] << "\"" 
                            << " -> " 
                            << "\"" << marker_ids[ii] << "\"" 
                            << "[label =" << buffer << "]" << ';' << endl;
                }
            }
        }
        ss_out1 << "}";        ss_out1 << endl;    ss_out1 << endl;
        ss_out2 << "}";        ss_out2 << endl;    ss_out2 << endl;
        return make_pair(ss_out1.str(),ss_out2.str());
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void linkage_group_graph_representation::contract_into_scc() {
        scc_adj_list.clear();
        scc_adj_list.resize(num_scc);
        map<int,int> oid_2_cid;
        for (int i = 0; i < num_scc; i++) {
            for (vector<int>::iterator iter = sccs[i].begin();
                 iter != sccs[i].end();
                 ++iter) {
                 oid_2_cid[*iter] = i;
            }
        }
        // for every scc find its adjacent nodes
        for (int i = 0; i < num_scc; i++) {
            map<int, vector<pair<double,double> > > adjct_list;
            // for every vertex in the SCC
            for (vector<int>::iterator iter2 = sccs[i].begin();
                 iter2 != sccs[i].end();
                 ++iter2) {
                for (vector<adj_edge_type>::iterator iter3 = adjacency_list[*iter2].begin();
                     iter3 != adjacency_list[*iter2].end();
                     ++iter3) {
                    int to_id = iter3->to_id;
                    double weight = iter3->weight;
                    double confidence = iter3->confidence;
                    int to_cc_id = oid_2_cid[to_id];
                    if (adjct_list.find(to_cc_id) == adjct_list.end()) {
                        adjct_list[to_cc_id] = vector<pair<double,double> >();
                        adjct_list[to_cc_id].push_back(make_pair(weight, confidence));
                    } else {
                        adjct_list[to_cc_id].push_back(make_pair(weight, confidence));
                    }
                }
            }
            // Consolidate the edges
            for (map<int, vector<pair<double,double> > >::iterator iter4 = adjct_list.begin(); 
                 iter4 != adjct_list.end();
                 ++iter4) {
                int to_cc_id = iter4->first;
                if (to_cc_id == i) continue;
                bool redundent = false;
                int from_cc_id = i;
                double weight = 0;
                double confidence = 0;
                for (vector<pair<double,double> >::iterator iter5 = (iter4->second).begin();
                     iter5 != (iter4->second).end();
                     ++iter5) {
                    weight = iter5->first * iter5->second + weight;
                    confidence = iter5->second + confidence;
                }
                weight = weight / confidence;
                scc_edge_type tmp_edge = {from_cc_id, to_cc_id, weight, confidence, redundent};
                scc_adj_list[from_cc_id].push_back(tmp_edge);
            }
        }
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void linkage_group_graph_representation::minimize_scc(){
        typedef boost::adjacency_list<vecS,vecS,directedS,no_property,property< edge_weight_t, int> > Graph;
        Graph sub_G(num_scc);
        for (int ii = 0; ii < num_scc; ii++) {
            for (vector<scc_edge_type>::iterator iter1 = scc_adj_list[ii].begin(); 
                 iter1 != scc_adj_list[ii].end(); 
                 ++iter1) {
                add_edge(ii, iter1->to_id, sub_G);
            }
        }
        
        property_map <Graph, edge_weight_t >::type w = get(edge_weight, sub_G);
        graph_traits<Graph>::edge_iterator e, e_end;
        for (boost::tie(e, e_end) = edges(sub_G); e != e_end; ++e) {
            w[*e] = 1;
        }
        
        /*Run johnson's all pairs shortest-path algorithm to figure out the transitive closure*/
        vector<vector<int> > D(num_scc, vector<int>(num_scc, -1));

        johnson_all_pairs_shortest_paths(sub_G, D);
        
        bool** d = new bool*[num_scc];
        for (int ii = 0; ii < num_scc; ii++) {
            d[ii] = new bool[num_scc];
        }

        for (int ii =0 ; ii < num_scc; ii++) {
            for (int jj = 0 ; jj < num_scc; jj++) {
                if (D[ii][jj] < 0) {
                    cout << "ERROR! the distance is not as expected" << endl;
                    assert(false); // crash the program
                }
                if (D[ii][jj] < num_scc +1) {
                    d[ii][jj] = true;
                } else {
                    d[ii][jj] = false;
                }
            }
        }
        
        for (int ii = 0 ; ii < num_scc; ii++) {
            d[ii][ii] = false;
        }
        
        vector<vector<int> > two_path(num_scc, vector<int>(num_scc, 0));
        for (int ii = 0 ; ii < num_scc; ii++) {
            for (int jj = 0; jj < num_scc; jj++) {
                for (int kk = 0; kk < num_scc; kk++) {
                    if (d[ii][kk] and (d[kk][jj])) {
                        two_path[ii][jj] = 1;
                        break;
                    }
                }
            }
        }
        
        /*delete the redundent edges*/
        num_redundent_scc_edge = 0 ; 
        for (int ii = 0; ii < num_scc; ii ++) {
            for (vector<scc_edge_type>::iterator iter1 = scc_adj_list[ii].begin(); 
                 iter1 != scc_adj_list[ii].end(); 
                 ++iter1) {
                if (two_path[ii][iter1->to_id] == 1) {
                    iter1->redundent = true;
                    num_redundent_scc_edge = num_redundent_scc_edge + 1;
                }
                if (iter1->from_id == iter1->to_id) {
                    iter1->redundent = true;
                }
            }
        }
        
        /*clean up*/
        for (int ii = 0 ; ii < num_scc ;ii++)
        {
            delete[] d[ii];
        }
        delete[] d;        
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    string linkage_group_graph_representation::output_dot_scc(LG_CLUSTER* lg_clster) {
        // To draw in the forward direction
        vector<SCC_Drawer*> drawers_forward(num_scc);
        for (int ii = 0; ii < num_scc; ii++) {
            vector<string> cur_scc;
            for (vector<int>::iterator iter1 = sccs[ii].begin(); iter1 != sccs[ii].end(); ++iter1) {
                cur_scc.push_back(marker_ids[*iter1]);
            }
            drawers_forward[ii] = lg_clster->GenDrawer(ii, true, cur_scc);
        }
        
        ostringstream ss_out(ios_base::out);
        ss_out << "digraph G {" << endl;
        
        // draw the nodes
        for (int ii = 0; ii < num_scc; ii++) {
            ss_out << drawers_forward[ii]->draw();
        }
        
        // draw the edges
        for (int ii = 0; ii < num_scc; ii++) {
            for (vector<scc_edge_type>::iterator iter2 = scc_adj_list[ii].begin();
                 iter2 != scc_adj_list[ii].end();
                 ++iter2) {
                if (iter2->redundent == false) {
                    int from_id = iter2->from_id;
                    int to_id = iter2->to_id;
                    double weight = iter2->weight;
                    char weight_buffer[20];
                    sprintf(weight_buffer, "%.2f", weight);
                    // for the forward direction
                    ss_out << "\"" << drawers_forward[from_id]->to() << "\"";
                    ss_out << " -> ";
                    ss_out << "\"" << drawers_forward[to_id]->from() << "\"";
                    ss_out << "[ label=\"" << weight_buffer << "\"]\n";
                    
                }
            }
        }        
        
        ss_out << "}" << endl << endl;

        
        for (int ii = 0; ii < num_scc; ii++) {
            delete drawers_forward[ii];
        }
        return ss_out.str();
    } 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    linkage_group_graph_representation::~linkage_group_graph_representation()
    {
        return;
    };
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
}
