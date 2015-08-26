#ifndef LARLITE_CLUSTERSHOWER_CXX
#define LARLITE_CLUSTERSHOWER_CXX

#include "ClusterShower.h"


namespace larlite {

  bool ClusterShower::initialize() {
   return true;
  }
  
  bool ClusterShower::analyze(storage_manager* storage) {
//Simple way to keep track for couts
std::cout<<"+++++++++++START++++++++++++++"<<std::endl;
Eventcounter++;
std::cout<<"Event number "<<Eventcounter<<std::endl;
//=================================================
//Define some things 
	// Need to clean this up
//=================================================
        auto geom = larutil::GeometryUtilities::GetME();
        auto tservice = larutil::TimeService::GetME();
	auto const& tpc_clock = tservice->TPCClock();
	double tick_offset = tservice->TriggerOffsetTPC() * tpc_clock.Frequency();// Kazu's fix
        unsigned int nplanes = larutil::Geometry::GetME()->Nplanes();
        // these things will be filled and used 
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$Define some variables$$$$$$$$$$$$$$$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

//=========================================
//========== Bring in info  ===============
//=========================================
	auto ev_track = storage->get_data<event_track>("cctrack");
	if(!ev_track){ print(msg::kERROR,__FUNCTION__,"No cclustertrack data product"), throw std::exception();}
	if(ev_track->empty()){ std::cout<<"No cctracks in this event"<<std::endl; return false;}

	event_hit* ev_hit = nullptr;
	auto const& track_hit_ass = storage->find_one_ass( ev_track->id(),
                                                       ev_hit,
                                                       ev_track->id().second );
	if(!ev_hit) std::cout<<"No hit association found..."<<std::endl;
	// Calo if we have? 
	//event_calorimetry* ev_calo = nullptr;
	//auto const& track_calo_ass = storage->find_one_ass( ev_track->id(),
         //                                               ev_calo,
         //                                               ev_track->id().second );
	//if(!ev_calo) std::cout<<"No calo association found..."<<std::endl;
	//

	//auto all_ev_hits = storage->get_data<event_hit>("cccluster");
	auto all_ev_hits = storage->get_data<event_hit>("gaushit");
	if(!all_ev_hits){ print(msg::kERROR,__FUNCTION__,"No cclusterhit data product"), throw std::exception();}
	if(all_ev_hits->empty()){ std::cout<<"No all cchits in this event"<<std::endl; return false;}

//------------------------------------------------------------
//%%%%------ 1. GET TRACKS AND LOOK AT DEDX OF START OF TRACK that has a good DEDX value for a photon.
//%%%%------ 2. Get the start point and guess direction. 
//%%%%------ 3. See if the projection gets filled with hits
//%%%%---------- If if does not then: 
//%%%%------------------------------ A. it might be going the wrong way   
//%%%%------------------------------ B. it might be a small shower   
//%%%%------------------------------ C. it might have a bad direction    
//%%%%------------------------------ C. it might be a track from proton,backscatter, ect...  
//%%%%------ 4. If the projection gets filled with hits
//%%%%------    Make the cone projection and see if there are any tracks inside of that projection 
//%%%%------    Add tracks to usedtracklets
//%%%%------   Store the tracklet and cone info somewhere 
//%%%%------ 5. ReCluster hits inside the cone and pass along tracklet info as axis 
	 

TVector3 BestTrackletPos;
TVector3 BestTrackletDir;
std::vector<int> UsedTracklets;
	// Make sure that we have atleast one tracklet
	for(size_t track_index = 0; track_index < ev_track->size(); track_index++) {
	//Check to see if the track is not already accounted for
		bool usedtrk = false;
		for(unsigned int check=0; check<UsedTracklets.size();check++) if(track_index == UsedTracklets[check]) usedtrk=true;
		if(usedtrk) continue;

			std::cout<<"Trajectory Points "<< ev_track->at(track_index).NumberTrajectoryPoints()<<std::endl;
			std::cout<<"Total Track Length "<< ev_track->at(track_index).Length()<<std::endl;
	//%%%%%--- 1Get the DEDX 
		// This should come from calo for now we will just try to calculate it here... 
		auto const& hit_index_v = track_hit_ass[track_index];
		std::cout<<" Size of hit vector "<<hit_index_v.size()<<std::endl;
		double dedxAvgCol = 0;// For now we just use colection... we can change this later and grab from calo
		double ColCounter = 0;
		for(unsigned int b=0; b<hit_index_v.size();b++){
			auto thehit = ev_hit->at( hit_index_v[b] ) ;
			unsigned int theplane = thehit.View();
			if(theplane==2){ ColCounter++;// I think 2 is the collection.... need to check this
			dedxAvgCol += fCaloAlg.dEdx_AMP(thehit.Integral(),(thehit.PeakTime()+ tick_offset)*geom->TimeToCm(), 3.0, theplane );
				}// if collection plane
			// What does pitch mean.... isn;t it just 3 cm
			}// for loop over associated hits
		//DEDX avg 
		dedxAvgCol/=ColCounter;
		std::cout<<"dEdx avg value for Collection :"<<dedxAvgCol<<std::endl;
		
		// If DEDX is good then use this... else continue 

	//%%%%%--- 2. Get the start of the track 
		BestTrackletPos = ev_track->at(track_index).Vertex();
		BestTrackletDir = ev_track->at(track_index).VertexDirection();

		for(size_t pt=0; pt< ev_track->at(track_index).NumberTrajectoryPoints(); pt++){
			 std::cout<<"\t Position at point ("<< pt << ")   : ["
				 <<ev_track->at(track_index).LocationAtPoint(pt)[0]<<","
				 <<ev_track->at(track_index).LocationAtPoint(pt)[1]<<" , "
				 <<ev_track->at(track_index).LocationAtPoint(pt)[2]<<" ]"<<std::endl;
			 std::cout<<"\t Direction at point ("<< pt << ")   : ["
				 <<ev_track->at(track_index).DirectionAtPoint(pt)[0]<<","
				 <<ev_track->at(track_index).DirectionAtPoint(pt)[1]<<" , "
				 <<ev_track->at(track_index).DirectionAtPoint(pt)[2]<<" ]"<<std::endl;
			}// Loop over all the points

		}// for loop over tracks
	
	
	

//%%%%------ 3. See if the projection gets filled with hits
	std::cout<<"&&&&CONE STUFF&&&&& "<<std::endl;
	auto coneintpc = fgeoconic.ConeInTPC(BestTrackletPos,BestTrackletDir,ConeLength,angle, smoothness);
	std::cout<<"Is the cone in the TPC? "<<coneintpc<<std::endl;
	// Even if the conse is not in the TPC we still can force it to fit
		// This can be cleaner but for now it's fine
	auto ConeEdge0 = fgeoconic.FitConicalFeatures(BestTrackletPos,BestTrackletDir,ConeLength,angle, 0,smoothness);
	auto ConeEdge1 = fgeoconic.FitConicalFeatures(BestTrackletPos,BestTrackletDir,ConeLength,angle, 1,smoothness);
	auto ConeEdge2 = fgeoconic.FitConicalFeatures(BestTrackletPos,BestTrackletDir,ConeLength,angle, 2,smoothness);
	
	std::cout<<std::endl;
		// This is useful if you just want to plot the polygon quickly
		for(int ep=0; ep<ConeEdge2.size(); ep++)
			std::cout<<ConeEdge2[ep].w<<" "<<ConeEdge2[ep].t+250.0<<std::endl;
	std::cout<<std::endl;
	// CHECK if this at all makes sense for a shower.
		// FOR NOW to just fill the code out.... lets just take the direction from the track. 


	
//%%%%------ 4. If the projection gets filled with hits and is deemed good 
		// Address the other tracklets that are in the volume
	
        for(size_t track_index = 0; track_index < ev_track->size(); track_index++) {
                auto testPos = ev_track->at(track_index).Vertex();
                std::vector<double> pos;
                pos.push_back(testPos.X());
                pos.push_back(testPos.Y());
                pos.push_back(testPos.Z());
		// just pic a plane to test
		auto test2d = geom->Get2DPointProjectionCM(pos,2);
		bool TrackInPoly = fgeoconic.TrackStartContain(test2d, ConeEdge2);
		if(TrackInPoly) UsedTracklets.push_back(track_index);	
		}

//%%%%------ 5. ReCluster hits inside the cone and pass along tracklet info as axis 

	// This can be done in the beginning. 
        //== Make the pxhit  vector by plane for now...
                std::vector<std::vector<larutil::PxHit>> PxHitsVect(nplanes);

		// Sort out the hits
		for(auto const& hit : *all_ev_hits){ 
                              ::larutil::PxHit h;
                              h.t = (hit.PeakTime() + tick_offset )* geom->TimeToCm() ;
                              h.w = hit.WireID().Wire     * geom->WireToCm();
                              h.charge = hit.Integral();
                              h.plane  = hit.View();
                                if( (int)hit.View() ==0) PxHitsVect[0].push_back(h);
                                if( (int)hit.View() ==1) PxHitsVect[1].push_back(h);
                                if( (int)hit.View() ==2) PxHitsVect[2].push_back(h);
                                  }

		std::vector<larutil::PxHit> contp0 = fgeoconic.PolyContain(PxHitsVect[0], ConeEdge0);
		std::vector<larutil::PxHit> contp1 = fgeoconic.PolyContain(PxHitsVect[1], ConeEdge1);
		std::vector<larutil::PxHit> contp2 = fgeoconic.PolyContain(PxHitsVect[2], ConeEdge2);
	
	// Needs to be a way to remove the hits from the PxHitsVect once they have been used in a cluster. 
	// Or we can assign then to each cone and make a choice later on before making the clusters.
	std::cout<<"&&&&HIT CLUSTER STUFF&&&&& "<<std::endl;
	// Lets see what this looks like 
	for( unsigned int a =0 ; a<contp2.size(); a++){
		std::cout<<contp2[a].w<<" "<<contp2[a].t+250.0<<std::endl;
		}
//------------------------------------------------------------



   // First of all create an output
//    auto Output_cluster = storage->get_data<event_cluster>("ClusterShower");









/*
        for(auto const& h : *hits) hitsvect.push_back(h);
                        for(auto const& hit : *hits){ 
                              ::larutil::PxHit h;
                              h.t = (hit.PeakTime() + tick_offset )* geom->TimeToCm() ;
                              h.w = hit.WireID().Wire     * geom->WireToCm();
                              h.charge = hit.Integral();
                              h.plane  = hit.View();
                                if( (int)hit.View() ==0) PxHitsVect[0].push_back(h);
                                if( (int)hit.View() ==1) PxHitsVect[1].push_back(h);
                                if( (int)hit.View() ==2) PxHitsVect[2].push_back(h);
                                  }
//  This uses the truth info to get the cone
        auto mcshower = storage->get_data<event_mcshower>("mcreco");
			for(auto const& mcs : *mcshower){
				auto SP = mcs.Start();
                                auto ShowerDetProf =  mcs.DetProfile();
				//StartConePos = SP.Position();
				//StartShowerPos = ShowerDetProf.Position();
					StartShowerPos = SP.Position();
					StartConePos = ShowerDetProf.Position();
				auto pos = ShowerDetProf.Position();
				//StartConeDir = SP.Momentum();
				//StartShowerDir = ShowerDetProf.Momentum();
					StartShowerDir = SP.Momentum();
					StartConeDir = ShowerDetProf.Momentum();
				auto dir = ShowerDetProf.Momentum();
				DetEn = ShowerDetProf.E();
				mcEn = SP.E();
				mctcontain = DetEn/mcEn;
						}//mcshower 
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-------------Bring in info-----------------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================


//=--=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
		std::vector<double>  ChargeinCone(nplanes);
		//coneintpc = fgeoconic.ConeInTPC(StartConePos,StartConeDir,ConeLength,angle, smoothness);
		bool showerintpc =  false;
		try{
		showerintpc = fgeoconic.ConeInTPC(StartShowerPos,StartShowerDir,ConeLength,angle, smoothness);
		showerintpc = fgeoconic.ConeInTPC(StartShowerPos,StartShowerDir,2,angle, smoothness);
		}catch(const ::larutil::LArUtilException& e){
		return false;
		showerintpc = false;
                }
		std::vector<double> ratio(nplanes);

	//if(coneintpc && energy>100 && energy<205)
	//if(showerintpc && DetEn/mcEn > 0.95 &&mcEn>200 &&mcEn<250)
	//if(showerintpc && DetEn/mcEn > 0.9 &&mcEn>200 &&mcEn<250 && StartConeDir.Pz()/StartConeDir.E()>0.8 &&StartConeDir.Px()/StartConeDir.E()<0.01)
	if(showerintpc && DetEn/mcEn > 0.9 &&mcEn>200 &&mcEn<255 && StartConeDir.Pz()/StartConeDir.E()>0.8&&StartConeDir.Px()/StartConeDir.E()>0.1 )
	//if(showerintpc &&mcEn<400)
	//asdf
	{
	//std::cout<<"pz" <<StartConeDir.Pz()/StartConeDir.E()<<std::endl;
		std::vector<double> DChargeDL(nplanes,0.0);
		for(unsigned int plane=0; plane<nplanes; plane++)
		{
			double DL = 2;
			for(double a=2.0; a<ConeLength ; a+=DL)
			{
			//auto polyproj = fgeoconic.ConicalFeatures(StartConePos, StartConeDir,a, angle, plane ,smoothness);
			std::vector<larutil::PxHit> contp;
			try{
			auto polyproj = fgeoconic.ConicalFeatures(StartShowerPos, StartShowerDir,a, angle, plane ,smoothness);
			contp = fgeoconic.PolyContain(PxHitsVect[plane], polyproj);
			}catch(const ::larutil::LArUtilException& e){
			continue;}
			// change the dir
				double newchargecont=0.0;
				for(auto const h: contp) newchargecont+=h.charge;
				double dCdX = newchargecont/a;
				if(dCdX<100 && a<14 ) continue;	
			}
		}
	}//
	
*/
std::cout<<"+++++++++++END++++++++++++++++"<<std::endl;
    return true;
  }

  bool ClusterShower::finalize() {

	if(_fout)
	_fout->cd();
    return true;
  }

}
#endif
