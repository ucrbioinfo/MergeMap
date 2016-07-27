/*
 *  fragment.h
 *  consensusmap
 *
 *  Created by yonghui on 10/15/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */
 
#ifndef FRAGMENT_HEADER
#define FRAGMENT_HEADER

#include <vector>
#include <string>
#include <cassert>
#include <map>
#include "shortest_path.h"

using namespace std;

namespace consensus_map{
    struct node_id_strct{
        int frag_id;
        int level_id;
        int mid;
    };
    
    struct ID_cmp{
        bool operator()(const node_id_strct & id1, const node_id_strct & id2) const {
            if (id1.frag_id < id2.frag_id) return true;
            else if (id1.frag_id > id2.frag_id) return false;
            else if (id1.level_id < id2.level_id) return true;
            else if (id1.level_id > id2.level_id) return false;
            else if (id1.mid < id2.mid) return true;
            else return false;
        }
    };
    
    class Fragment{
        public:
            // Some simple accessors. They are all inline functions
            int get_num_bins() const{
                return markers_in_bins_.size();
            };
            
            const vector<vector<string> > & get_bins() const {
                return markers_in_bins_;
            };
            
            const vector<double> & get_distances() const {
                return distance_btw_adjacent_pairs_;
            }
            
            double get_confidence() const {
                return confidence_;
            }
            
            string get_description() const {
                return description_;
            }
            
            int get_num_mrks() const {
                int num_mrks = 0;
                for (int ii = 0; ii < markers_in_bins_.size(); ii++) {
                    num_mrks = num_mrks + markers_in_bins_[ii].size();
                }
                return num_mrks;
            }
            
            bool get_prob_set() const {
                return prob_set_;
            }
            
            const vector<vector<double> > & get_prob() const {
                return prob_del_;
            }
            
            void set_prob(vector<vector<double> > & _in_prob) {
                prob_del_ = _in_prob;
                prob_set_ = true;
            }

            bool get_del_set() const {
                return del_set_;
            }
            
            const vector<vector<bool> > & get_heuristic_del() const {
                return heuristic_del_;
            }
            
            const vector<vector<bool> > & get_del() const {
                return dels_;
            }
            
            void set_del(vector<vector<bool> > & _in_del) {
                dels_ = _in_del;
                del_set_ = true;
            }
            
            void h_delete(int bin_id, int marker_id){
                assert(bin_id < markers_in_bins_.size());
                assert(bin_id >= 0);
                assert(marker_id < markers_in_bins_[bin_id].size());
                assert(marker_id >= 0);
                assert(heuristic_del_[bin_id][marker_id] == false);
                heuristic_del_[bin_id][marker_id] = true;
            };
            
            void FlipOrientation();
            
            // Define the constructor, copy constructor and assignment operator
            Fragment(const vector<vector<string> >& _markers_in_bins, 
                     const vector<double>& _distance_btw_adjacent_pairs,
                     double _confidence,
                     string _description){
                markers_in_bins_ = _markers_in_bins;
                distance_btw_adjacent_pairs_ = _distance_btw_adjacent_pairs;
                confidence_ = _confidence;
                description_ = _description;
                prob_set_ = false;
                del_set_ = false;
                // initialize the heuristic_del_ object
                heuristic_del_.resize(markers_in_bins_.size());
                for (int ii = 0; ii < markers_in_bins_.size(); ii++) {
                    heuristic_del_[ii].resize(markers_in_bins_[ii].size(), false);
                }
            };
            
            Fragment(const Fragment& _frag){
                markers_in_bins_ = _frag.get_bins();
                distance_btw_adjacent_pairs_ = _frag.get_distances(); 
                confidence_ = _frag.get_confidence();
                description_ = _frag.get_description();
                prob_set_ = _frag.get_prob_set();
                prob_del_ = _frag.get_prob();
                del_set_ = _frag.get_del_set();
                dels_ = _frag.get_del();
                heuristic_del_ = _frag.get_heuristic_del();
            };
            
            Fragment& operator=(const Fragment& _frag){
                markers_in_bins_ = _frag.get_bins();
                distance_btw_adjacent_pairs_ = _frag.get_distances();
                confidence_ = _frag.get_confidence();
                description_ = _frag.get_description();
                prob_set_ = _frag.get_prob_set();
                prob_del_ = _frag.get_prob();
                del_set_ = _frag.get_del_set();
                dels_ = _frag.get_del();
                heuristic_del_ = _frag.get_heuristic_del();
                return *this; 
            };
            
            // Use the default destructor
            ~Fragment(){ 
                return;
            }
            
            vector<string> return_deleted(){
                vector<string> to_return;
                for (int ii = 0; ii < markers_in_bins_.size(); ii++) {
                    for (int jj = 0; jj < markers_in_bins_[ii].size(); jj++) {
                        if (dels_[ii][jj]) {
                            to_return.push_back(markers_in_bins_[ii][jj]);
                        }
                        if (heuristic_del_[ii][jj]) {
                            to_return.push_back(markers_in_bins_[ii][jj]);
                        }
                    }
                }
                return to_return;
            };
            
            void dump();
            
        private:
            // The individual map in the linear order of bins, where each bin is a set of markers 
            // The distances between adjacent bins are stored in the distance_btw_adjacent_pairs_
            // The size of the distance_btw_adjacent_pairs_ vector is 1 minus the size of the number of bins
            // The distance vector is empty if the number of bins is equal to 1
            // The distance is in CMs
            vector<vector<string> > markers_in_bins_; 
            vector<double> distance_btw_adjacent_pairs_; 
            
            // Confidence in the accuracy of the fragment.
            // The purpose of having the confidence here is that when assembling fragments of conflicting order,
            // those markers from less confident fragment will get deleted more likely
            double confidence_;
            string description_;
            
            // prob_del stores the probability to delete a particular node in order to solve conflicts 
            // among related fragments
            bool prob_set_;
            vector<vector<double> > prob_del_;
            
            bool del_set_;
            vector<vector<bool> > dels_;
            
            vector<vector<bool> > heuristic_del_;
        
    };
    
    // When building consensus maps, individual map is chopped into pieces called fragments.
    // When two fragments share markers, they are related.
    // The following class is intended to encapsulate the operations needed to assemble related fragments 
    // into a consensus order.
    class FragClstr{
        public:
            // several accessors
            const vector<Fragment>& GetFrags() const {
                return frags_;
            }
            const vector<int>& GetFragIds() const {
                return frag_ids_;
            }
            
            // Constructrs 
            FragClstr() {
                return;
            }
            void Initialize(const vector<Fragment>& _frags, const vector<int>& _frag_ids) {
                frags_ = _frags;
                frag_ids_ = _frag_ids;
            };
            
            FragClstr& operator=(const FragClstr& rhs) {
                frags_ = rhs.GetFrags();
                frag_ids_ = rhs.GetFragIds();
                return (*this);
            }
            
            FragClstr(const FragClstr& rhs) {
                (*this) = rhs;
            };
            
            int total_num_mrks(){
                int total_mrks = 0;
                for ( int ii = 0; ii < frags_.size(); ii++) {
                    total_mrks = total_mrks + frags_[ii].get_num_mrks();
                }
                return total_mrks;
            }
            
            void FlipOrientation() {
                for (int ii =0; ii < frags_.size(); ii++) {
                    frags_[ii].FlipOrientation();
                }
            }

            // The following function will execute the lagrangian optimization algorithm to solve the LP
            void solveLP(double epsilon);
            
            vector<pair<string,int> > return_conflicts();
            
            void h_delete(int frg_id, int bin_id, int marker_id){
                assert(frg_id < frags_.size());
                assert(frg_id >= 0);
                frags_[frg_id].h_delete(bin_id, marker_id);
            };
            
            void dump();
            
        private:
            // The following function will populate the graph_lag object
            CG* GenGraphLag(map<int, node_id_strct>& node_id_map);
            

            // The fragments to be assembled are stored a vector
            // The individual fragments are assumed to be in the right order already
            vector<Fragment> frags_;
            vector<int> frag_ids_;
            
    };
    
    // In the process of building consensus map, each individual maps are to be chopped into 
    // pieces called fragments. 
    // The following class stores the sequence of fragments for one single population.
    // The distance between adjacent fragments are stored in the distances_ vector.
    // The size of the distance_ vector is 1 minus the size of the frags_ vector.
    class FragSeq{
        public:
            FragSeq(const vector<Fragment>& _frags, const vector<double>& _distances){
                frags_ = _frags;
                distances_ = _distances;
            }
        private:
            vector<Fragment> frags_;
            vector<double> distances_;
    };
};

#endif
