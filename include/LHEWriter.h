#include <fstream>
#include <vector>
#include <string>

#include "particles.h"
#include "TLorentzVector.h"

using namespace std;

class LHEWriter
{
    private:
        ofstream oF;

    public:
        LHEWriter(string fName, double sqrtS, bool isweight);
        ~LHEWriter();

        int writeEvent(vector<particle> outParts, double Q, double weight = 1 );
        int close();

};
