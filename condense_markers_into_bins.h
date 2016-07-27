/*
 *  condense_markers_into_bins.h
 *  consensusmap
 *
 *  Created by yonghui on 5/12/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONDENSE_MARKERS_INTO_BINS
#define CONDENSE_MARKERS_INTO_BINS
 
#include <vector>
#include <set>
#include <string>
#include <utility>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <boost/regex.hpp>
#include <boost/graph/strong_components.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/johnson_all_pairs_shortest.hpp>
#include "constants.h"
#include <map>
#include <utility>
#include "single_population.h"

using namespace std;
using namespace boost;
namespace consensus_map
{
    class linkage_group_bin;
    class linkage_group_DAG;
    // _in_p_to_array_lgs is a reference to an array of pointers to linkage_group_bin objects
    // _out_p_to_array_lgs is the output lgs where the markers are condensed into supernodes
    // The individual elements of the array is initialized inside the function, 
    // and should be deallocated by the caller of this function. 
    void consolidate_markers_into_bins(vector<const linkage_group_bin*>& _in_p_to_array_lgs, 
                                       int _in_number_of_lgs, 
                                       vector<linkage_group_bin*>& _out_p_to_array_lgs, 
                                       const linkage_group_DAG & lg_dag);
};

#endif
