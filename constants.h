/*
 *  constants.h
 *  consensusmap
 *
 *  Created by yonghui on 5/16/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CONSTANTS_HEADER
#define CONSTANTS_HEADER

#include <string>
using namespace std;

const string MRKS_IN_LGS_RE = "The size of the linkage groups are:\\s+((?:\\d+\\s+)+)";
const string BINS_IN_LGS_RE = "The number of bins in each linkage group:\\s+((?:\\d+\\s+)+)";
//const string LGS_RE = "=+STATs for linkage group\\s*(\\d+\\.\\d+|\\d+)=+\\s*lowerbound:(\\d+\\.\\d+|\\d+)\\s*upperbound:\\s*(\\d+\\.\\d+|\\d+)\\s*cost after initialization:(?:\\d+\\.\\d+|\\d+)\\s*((?:\\s*\\w+\\s*\\d+\\s*)+)\\s*distances between adjacent pairs\\s*((?:(?:\\d+|\\d+\\.\\d+)\\s+)*)";
const string LGS_RE = "START OF LG\\s*"
                      "=+STATs for linkage group\\s*(\\d+)=+\\s*"
                      "lowerbound:(\\d+\\.\\d+|\\d+)\\s*"
                      "upperbound:\\s*(\\d+\\.\\d+|\\d+)\\s*"
                      "cost after initialization:(?:\\d+\\.\\d+|\\d+)\\s*"
                      "((?:\\s*\\w+\\s*\\d+\\s*)+)\\s*"
                      "distances between adjacent pairs\\s*((?:(?:\\d+|\\d+\\.\\d+)\\s+)*)";

const int MAXIMUM_NUM_OF_LGS_PER_MAP = 1000;
const double EPSILON = 0.000001;
const int FILE_NAME_MAX_LENGTH = 50;
const int DEFAULT_EDGE_CONFIDENCE = 1;
const double MAXIMUM_DISTANCE_IN_CM = 1000000;
const int MaxMrksPLine = 5;
const int CFG_LINE_LENGTH_MAX = 500;
const double LP_epsilon = 0.05;
const string start_lg = ";BEGINOFGROUP";
const string end_lg = ";ENDOFGROUP"; 

// define some constants for drawing dags
const double kHue = 1.0;
const double kSaturation = 0.0;
const double kBrightness = 0.8;
const int kLargestSCCLP = 50;


const int klgMinCommon = 5;
#endif
