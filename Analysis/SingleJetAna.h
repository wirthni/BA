#ifndef METHODS_H
#define METHODS_H
#include <Math/Vector4D.h>
#include <unordered_map>
#include <vector>
#include "functions.h"
#include "StJetTrackEvent.h"
#include "StJetTrackEventLinkDef.h"
#include "TClassRef.h"
#include <algorithm>
#include <cstdlib>
#include "fastjet/config.h"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "ReadTree.h"
#endif


using namespace fastjet;
using namespace std;

vector<PseudoJet> FindJets(const vector<PseudoJet> vec_particles);

class MyUserInfo : public fastjet::PseudoJet::UserInfoBase {
public:
    MyUserInfo(int NoOfParticles) 
        : m_NoOfParticles(NoOfParticles){}

    int getNoOfParticles() const { return m_NoOfParticles; }

private:
    int m_NoOfParticles;
};
