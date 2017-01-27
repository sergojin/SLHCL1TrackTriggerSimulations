#ifndef AMSimulationIO_CSVFileReader_h_
#define AMSimulationIO_CSVFileReader_h_

#include <map>
#include <vector>

#include "TString.h"

class CSVFileReader {
  public:
    CSVFileReader(int verbose=1);
    ~CSVFileReader();

    void getTriggerTowerMap(TString src, unsigned tt,
                            std::map<unsigned, std::vector<unsigned> >& ttmap,
                            std::map<unsigned, std::vector<unsigned> >& ttrmap);

  protected:
    const int verbose_;
};

#endif

