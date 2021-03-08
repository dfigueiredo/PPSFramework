#include "GenParticlesEvent.h"

// Retrieve GenParticles
std::vector<const reco::GenParticle*> GenParticlesEvent::GetGenParticles(){
	return genparticleslist;
}

// Class Definition, loop over all genparticles events and fill vector.
GenParticlesEvent::GenParticlesEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, edm::EDGetTokenT<std::vector< reco::GenParticle>>& GenPartToken_){

	genparticleslist.clear();

	try{
		// Get the GenParticle collection from the event
		edm::Handle<std::vector<reco::GenParticle>> GenPartColl; 
		iEvent.getByToken( GenPartToken_, GenPartColl );
		for (unsigned int i = 0; i < GenPartColl->size(); ++i ) {
			const reco::GenParticle* genpart = &((*GenPartColl)[i]);
			genparticleslist.push_back(genpart);
		}
		LogDebug( "GenParticles Event" ) << "Passed Loop on genparticle";
	}catch(...){
		LogDebug( "GenParticles Event" ) << "collection not present in the sample";
	}

}

