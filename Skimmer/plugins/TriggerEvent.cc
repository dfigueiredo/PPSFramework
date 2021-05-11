// -*- C++ -*-
//
// Package:    CTPPSAnalysisCode/TriggerEvent
// Class:      TriggerEvent
// 
/**\class TriggerEvent TriggerEvent.cc CTPPSAnalysisCode/TriggerEvent/plugins/TriggerEvent.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Diego Figueiredo
//         Created:  Fri, 05 Oct 2018 09:08:27 GMT
//
//

#include "TriggerEvent.h"

// Retrieve Trigger
std::vector<int> TriggerEvent::GetTrigger(){
  return triggerlist;
}

// L1 prescale
std::vector<int> TriggerEvent::GetPrescalesL1(){
  return  prescalesL1;
}

// HLT prescale
std::vector<int> TriggerEvent::GetPrescalesHLT(){
  return  prescalesHLT;
}

// Class Definition, loop over all trigger events and fill vector.
TriggerEvent::TriggerEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, HLTPrescaleProvider &hltPrescaleProvider, edm::EDGetTokenT<edm::TriggerResults>& triggerResultsToken_, std::vector<std::string>& triggersList_, bool enablePrescales_){

  bool debug = false;
  triggerlist.clear();

  try{
    edm::Handle<edm::TriggerResults> hltResults;
    iEvent.getByToken( triggerResultsToken_, hltResults);

    if( hltResults.isValid() ){

      unsigned int nSize = hltResults->size();
      const edm::TriggerNames& triggerNames = iEvent.triggerNames(*hltResults);

      size_t idxpath = 0;
      std::vector<std::string>::const_iterator hltpath = triggersList_.begin();
      std::vector<std::string>::const_iterator hltpaths_end = triggersList_.end();

      for(; hltpath != hltpaths_end; ++hltpath,++idxpath){

	std::string resolvedPathName;
	if( edm::is_glob( *hltpath ) ){
	  std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(triggerNames.triggerNames(), *hltpath);
	  if( matches.empty() ){
	    LogDebug("Configuration") << "Could not find trigger " << *hltpath << " in the path list.\n";
	  } 
	  else if( matches.size() > 1)
	    throw cms::Exception("Configuration") << "HLT path type " << *hltpath << " not unique\n";
	  else resolvedPathName = *(matches[0]);
	} else{
	  resolvedPathName = *hltpath;
	} 

	unsigned int idx_HLT = triggerNames.triggerIndex(resolvedPathName);

	if (idx_HLT < nSize){
	  int accept_HLT = ( hltResults->wasrun(idx_HLT) && hltResults->accept(idx_HLT) ) ? 1 : 0;
	  triggerlist.push_back(accept_HLT);
	  if(enablePrescales_){
	    const std::pair<int,int> prescales(hltPrescaleProvider.prescaleValues(iEvent,iSetup,resolvedPathName));
	    prescalesL1.push_back(prescales.first);
	    prescalesHLT.push_back(prescales.second);
	    if(debug){
	      std::cout << "--> Prescales L1T, HLT: " << prescales.first << ", " << prescales.second << std::endl;
	      const std::pair<std::vector<std::pair<std::string,int> >,int> prescalesInDetail(hltPrescaleProvider.prescaleValuesInDetail(iEvent,iSetup,resolvedPathName));
	      std::ostringstream message;
	      for (unsigned int i=0; i<prescalesInDetail.first.size(); ++i) {
		message << " " << i << ":" << prescalesInDetail.first[i].first << "/" << prescalesInDetail.first[i].second;
	      }
	      std::cout << "--> Detailed Info: " << resolvedPathName << " [" << idx_HLT << "] "
		<< std::endl
		<< "Prescales L1T: " << prescalesInDetail.first.size() <<  message.str()
		<< std::endl
		<< ", Prescale HLT: " << prescalesInDetail.second
		<< std::endl;
	    }
	  }else{
	    prescalesL1.push_back(-1);
	    prescalesHLT.push_back(-1);
	  }	  
	}else{
	  triggerlist.push_back(-1);
	  prescalesL1.push_back(-1);
	  prescalesHLT.push_back(-1);
	} 
      } 
    }else{
      edm::LogWarning("Trigger Event") << "trigger is not valid.\n";
    } 
  }catch(...){
    triggerlist.push_back(-1);
    prescalesL1.push_back(-1);
    prescalesHLT.push_back(-1);
  }

}
