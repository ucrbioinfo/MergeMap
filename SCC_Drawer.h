/*
 *  SCC_Drawer.h
 *  consensusmap
 *
 *  Created by yonghui on 10/18/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

 
#ifndef SCC_DRAWER_HEADER
#define SCC_DRAWER_HEADER

#include <string>
#include "ConsensusLG.h"

using namespace std;
namespace consensus_map{
    class SCC_Drawer {
        public:
            virtual string from() = 0;
            virtual string to() = 0;
            virtual string draw() = 0;
            string BreakNode(string in_str);
        private:
    };
    
    class Sigleton_Drawer: public SCC_Drawer{
        public:
            Sigleton_Drawer(string mrk_id, int scc_id) {
                mrk_id_ = mrk_id;
                scc_id_ = scc_id;
            };
            virtual string from() {
                char buffer[100];
                sprintf(buffer, "%d", scc_id_);
                return string(buffer);
            };
            virtual string to() {
                char buffer[100];
                sprintf(buffer, "%d", scc_id_);
                return string(buffer);
            };
            virtual string draw() {
                char buffer[500];
                sprintf(buffer, "\"%d\"[label=\"%s\"]\n", scc_id_, BreakNode(mrk_id_).c_str());
                return string(buffer);
            };
        private:
            string mrk_id_;
            int scc_id_;
    };
    
    class Multi_Drawer: public SCC_Drawer{
        public:
            Multi_Drawer(const FragClstr& frags, int scc_id) {
                frags_ = frags;
                scc_id_ = scc_id;
            };
            virtual string from(){
                char buffer[100];
                sprintf(buffer, "%d_1", scc_id_);
                return string(buffer);                
            };
            virtual string to(){
                char buffer[100];
                sprintf(buffer, "%d_2", scc_id_);
                return string(buffer);                
            };
            virtual string draw();
        private:
            string draw_fragment(int frag_id);
            FragClstr frags_;
            int scc_id_;
    };
}
#endif

