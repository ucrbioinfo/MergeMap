/*
 *  SCC_Drawer.cpp
 *  consensusmap
 *
 *  Created by yonghui on 10/18/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#include "SCC_Drawer.h"
#include "fragment.h"
#include "constants.h"
#include <iostream>
#include <sstream>

namespace consensus_map {
    string SCC_Drawer::BreakNode(string in_str) {
        // insert endl at every MaxMrksPLine comma
        int length = in_str.length();
        vector<char> new_node;
        int comma_id = 0;
        for (int ii = 0; ii < length; ii++) {
            if (in_str.at(ii) == ',') {
                comma_id = comma_id + 1;
                if (comma_id % MaxMrksPLine == 0) {
                    new_node.push_back('\\');
                    new_node.push_back('n');
                } else {
                    new_node.push_back(in_str.at(ii));
                }
            } else {
                new_node.push_back(in_str.at(ii));
            }
        }
        char buffer[new_node.size() + 1];
        buffer[new_node.size()] = 0;
        for (int ii = 0; ii < new_node.size(); ii++) {
            buffer[ii] = new_node[ii];
        }
        return string(buffer);
    }
    string Multi_Drawer::draw_fragment(int frag_id){
        const vector<Fragment>& frags = frags_.GetFrags();
        const Fragment& crt_fragment = frags[frag_id];
        const vector<vector<string> >& mrks_bins = crt_fragment.get_bins();
        const vector<vector<double> >& del_probs = crt_fragment.get_prob();
        const vector<vector<bool> >& dels = crt_fragment.get_del();
        const vector<vector<bool> >& hdels = crt_fragment.get_heuristic_del();
        
        if (mrks_bins.size() <= 1) return ""; // don't draw fragment with only one level
        const vector<double>& distances = crt_fragment.get_distances();
        ostringstream ss_out(ios_base::out);
        char fid_buffer[20];
        sprintf(fid_buffer, "%d", frag_id);
        char cluster_id[20];
        sprintf(cluster_id, "%d_%d", scc_id_, frag_id);
        char label[100];
        sprintf(label, 
                "%s (%.2f)", 
                (crt_fragment.get_description()).c_str(), 
                crt_fragment.get_confidence());
        char default_color[100];
        sprintf(default_color, "\"%.3f %.3f %.3f\"", kHue, kSaturation, kBrightness);
        ss_out <<  "subgraph cluster_" << cluster_id <<"{\n" 
               <<  "  style=filled;\n"
               <<  "  label=\"" << label << "\";\n"
               <<  "  color=" << default_color << ";\n" ;
        // draw the nodes
        for (int ii = 0; ii < mrks_bins.size(); ii++) {
            for (int jj = 0; jj < mrks_bins[ii].size(); jj++) {
                string nid = "\"" + mrks_bins[ii][jj] + fid_buffer + "\"";
                char label_buf[mrks_bins[ii][jj].length() + 100];
                sprintf(label_buf, "%s\\n(%.2f)", mrks_bins[ii][jj].c_str(), del_probs[ii][jj]);
                double saturation = del_probs[ii][jj];
                if (hdels[ii][jj]) {
                    saturation = 1.0;
                }
                char fill_color[100];
                sprintf(fill_color, ", fillcolor=\"%.3f %.3f %.3f\"", kHue, saturation, kBrightness);
                string shape;
                if (dels[ii][jj]) {
                    shape = ", shape=Mdiamond";
                } else if (hdels[ii][jj]){
                    shape = ", shape=tripleoctagon";
                } else {
                    shape = "";
                }
                
                ss_out << nid << "[label=\"" << BreakNode(label_buf) << "\"" 
                       << shape << fill_color 
                       << ", style=filled" << "];\n";
            }
        }
        for (int ii = 1; ii < mrks_bins.size(); ii++) {
            for (vector<string>::const_iterator iter1 = mrks_bins[ii-1].begin();
                 iter1 != mrks_bins[ii-1].end();
                 ++iter1) {
                double weight = distances[ii-1];
                char weight_buffer[50];
                sprintf(weight_buffer, "%.2f", weight);
                for (vector<string>::const_iterator iter2 = mrks_bins[ii].begin();
                     iter2 != mrks_bins[ii].end();
                     ++iter2) {
                    string from_id = "\"" + *iter1 + fid_buffer + "\"";
                    string to_id = "\"" + *iter2 + fid_buffer + "\"";
                    ss_out << from_id << " -> " << to_id << "[label=\"" << weight_buffer << "\"];\n"; 
                }
            }
        }
        // draw the edges
        ss_out << "}\n";
        return ss_out.str();
    };

    string Multi_Drawer::draw() {
        const vector<Fragment>& frags = frags_.GetFrags();
        ostringstream ss_out(ios_base::out);
        char cluster_id[20];
        sprintf(cluster_id, "%d", scc_id_);
        ss_out <<  "subgraph cluster_" << cluster_id << "{\n" 
               <<  "  color=blue;\n" ;

        // draw two special nodes
        string header = "\"" + from() + "\"";
        string tail = "\"" + to() + "\"";
        ss_out << header << "[label=\"\", shape=ellipse, width=0.1, height=0.1];\n";
        ss_out << tail << "[label=\"\", shape=ellipse, width=0.1, height=0.1];\n";

        // draw each individual fragment
        for (int ii = 0; ii < frags.size(); ii++) {
            ss_out << draw_fragment(ii);
        }
            
        // draw edges from header and tail to the header and tail of each fragment
        for (int ii = 0; ii < frags.size(); ii++) {
            const Fragment& crt_fragment = frags[ii];
            const vector<vector<string> >& mrks_bins = crt_fragment.get_bins();
            if (mrks_bins.size() <= 1) continue; // don't draw fragment with one level
            char fid_buffer[20];
            sprintf(fid_buffer, "%d", ii);
            for (vector<string>::const_iterator iter1 = mrks_bins[0].begin();
                 iter1 != mrks_bins[0].end();
                 ++iter1) {
                string to_id = "\"" + *iter1 + fid_buffer + "\"";
                ss_out << header << " -> " << to_id << ";" << endl;
            }
            int crt_num_bins = mrks_bins.size();
            for (vector<string>::const_iterator iter1 = mrks_bins[crt_num_bins - 1].begin();
                 iter1 != mrks_bins[crt_num_bins - 1].end();
                 ++iter1) {
                 string from_id = "\"" + *iter1 + fid_buffer +  "\"";
                 ss_out << from_id << " -> " << tail << ";" << endl;
            }
        }
        ss_out << "}\n";
        return ss_out.str();
    };
}
