#include "SLHCL1TrackTriggerSimulations/AMSimulationIO/interface/CSVFileReader.h"

#include <cassert>
#include <stdexcept>
#include <fstream>
#include <sstream>


CSVFileReader::CSVFileReader(int verbose)
: verbose_(verbose) {}

CSVFileReader::~CSVFileReader() {}


void CSVFileReader::getTriggerTowerMap(TString src, unsigned tt,
                                       std::map<unsigned, std::vector<unsigned> >& ttmap,
                                       std::map<unsigned, std::vector<unsigned> >& ttrmap) {

    if (!src.EndsWith(".csv")) {
        TString msg = "Input filename must be .csv";
        throw std::invalid_argument(msg.Data());
    }

    //if (verbose_)  std::cout << Info() << "Opening " << src << std::endl;
    std::vector<unsigned> values;

    std::string line, line2;
    std::ifstream ifs(src.Data());

    unsigned i = 0;
    while (std::getline(ifs, line)) {
        std::istringstream iss(line);

        if (i == tt + 1) {  // skip the first line
            unsigned j = 0;
            while (std::getline(iss, line2, ',')) {  // split by commas
                if (j == 0) {
                    assert(int(tt / 8) + 1 == std::stoi(line2));
                } else if (j == 1) {
                    assert(int(tt % 8) + 1 == std::stoi(line2));
                } else {
                    values.push_back(std::stoi(line2));
                }

                ++j;
            }
        }
        ++i;
    }

    if (values.empty()) {
        TString msg = "Failed to read any trigger tower";
        throw std::invalid_argument(msg.Data());
    }

    // Fill ttmap
    ttmap.clear();
    ttmap.insert(std::make_pair(tt, values));

    // Fill ttrmap
    ttrmap.clear();
    std::vector<unsigned>::const_iterator it;
    for (it = values.begin(); it != values.end(); ++it) {
        ttrmap[*it].push_back(tt);
    }
}
