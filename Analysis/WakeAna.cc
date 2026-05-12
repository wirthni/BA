#include "./StJetTrackEvent.h"
#include "./StJetTrackEventLinkDef.h"
#include "./functions.h"
#include <algorithm>
using namespace std;
#if false
#include "fastjet/config.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Subtractor.hh"
#include "StPhysicalHelixD.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#endif

const char JETTRACK_EVENT_TREE[]   = "JetTrackEvent";
const char JETTRACK_EVENT_BRANCH[] = "Events";

