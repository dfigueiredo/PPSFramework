#include "GenProtonsEvent.h"

// Retrieve GenProtons
std::vector<const reco::GenParticle*> GenProtonsEvent::GetGenProtons(){
	return genparticleslist;
}

// Class Definition, loop over all genparticles events and fill vector.
GenProtonsEvent::GenProtonsEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<std::vector< reco::GenParticle>>& GenPartToken_){

	genparticleslist.clear();

	try{
		// Get the GenParticle collection from the event
		edm::Handle<std::vector<reco::GenParticle>> GenPartColl; 
		iEvent.getByToken( GenPartToken_, GenPartColl );
		for (unsigned int i = 0; i < GenPartColl->size(); ++i ) {
			const reco::GenParticle* genpart = &((*GenPartColl)[i]);
                        if(genpart->pdgId() == 2212 && genpart->status()==1) genparticleslist.push_back(genpart);
		}
		LogDebug( "GenProtons Event" ) << "Passed Loop on gen protons";
	}catch(...){
		LogDebug( "GenProtons Event" ) << "collection not present in the sample";
	}

}

