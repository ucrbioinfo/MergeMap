#include <iostream>
#include "single_population.h"
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "constants.h"
#include "map_integration.h"

using namespace std;
using namespace consensus_map;

struct MapConfig {
    string ppl_name;
    double confidence;
    string file_path;
};

int main (int argc, char * const argv[]) 
{

    char buffer[CFG_LINE_LENGTH_MAX];
    vector<MapConfig> maps_to_read;
    
    ifstream cfg_file(argv[1]);
    
    cfg_file.getline(buffer, CFG_LINE_LENGTH_MAX);
    
    // parse the config file into a list of MapConfig
    // read one line from the input file
    while (not cfg_file.fail()) {
        if (not (buffer[0] == '#')) {
            istringstream crt_map(buffer);
            string ppl_name;
            double confidence;
            string file_path;
            crt_map >> ppl_name;
            crt_map >> confidence;
            crt_map >> file_path;
            MapConfig tmp_cfg = {ppl_name, confidence, file_path};
            maps_to_read.push_back(tmp_cfg);
        }
        cfg_file.getline(buffer, CFG_LINE_LENGTH_MAX);        
    }
    cfg_file.close();
    
    printf("number of maps %d \n", (int)maps_to_read.size());
    for (int i = 0; i < maps_to_read.size(); i++) {
        printf("%s,\t%.3f,\t%s\n", 
               (maps_to_read[i]).ppl_name.c_str(), 
               (maps_to_read[i]).confidence, 
               (maps_to_read[i]).file_path.c_str());
    }

    int number_of_ppls = maps_to_read.size(); 
    individual_mapping_ppl* individual_maps = new individual_mapping_ppl[number_of_ppls];
    for (int ii = 0 ; ii < number_of_ppls; ii++) {
        individual_maps[ii].read_from_file(maps_to_read[ii].file_path, 
                                           maps_to_read[ii].confidence, 
                                           maps_to_read[ii].ppl_name);
    }

    // consensus_map object will take care of the individual_maps array
    // there is no need to call delete[] individual_mapping_ppl at the end of the main function
    map_integration consensus_map_barley(individual_maps, number_of_ppls);
    consensus_map_barley.dump();
    consensus_map_barley.generate_dot_files();
    consensus_map_barley.Remove_conflicts();
    consensus_map_barley.Gen_consensus_dot();
    consensus_map_barley.Linearize_Dag();
    consensus_map_barley.Gen_linear_graph();
    consensus_map_barley.Gen_map_chart();
    return 0;
}
