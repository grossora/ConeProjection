#ifndef LARLITE_CLUSTERSHOWER_CXX
#define LARLITE_CLUSTERSHOWER_CXX

#include "ClusterShower.h"


namespace larlite {

  bool ClusterShower::initialize() {
       // fmcEnGood = new TH1D("fmcEnGood","Good Showers ",200,0,1000);
   return true;
  }
  
  bool ClusterShower::analyze(storage_manager* storage) {
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
	//== LarLight vector of hits
		std::vector<larlite::hit> hitsvect;
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::pair<larutil::PxPoint,larutil::PxPoint>> AxisSEpt(nplanes);
	//== 2D: This had the vertex point and end point in two 2d 
		std::vector<std::vector<larutil::PxPoint>> ConePolygonProjection(nplanes);
	//== Slope and Cept of the cone axis
		std::vector<std::pair<double,double>> ConeAxisSlopeCept;
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<larutil::PxHit>> PxHitsVect(nplanes);
	//== RecoFit Based on weights 
		std::vector<std::pair<double,double>> recofitvec(nplanes);
	//== Truth Photon Start Position  
		TLorentzVector StartConePos;
		TLorentzVector StartShowerPos;
	//== Truth Photon Start Dir  
		TLorentzVector StartConeDir;
		TLorentzVector StartShowerDir;
        //== Make the pxhit  vector by plane for now... using hits
		std::vector<std::vector<larutil::PxHit>> ContainedPxHitsVect(nplanes);
        //== Make the pxhit  vector by plane for now...
		std::vector<std::vector<unsigned int>> ContainedHitsVect(nplanes);
        //== Detector energy...
		double DetEn = -999;
//=================================================
//$$$$$$$$$$$$$$$$$---END---$$$$$$$$$$$$$$$$$$$$$$$
//-----------Define some variables-----------------
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//=================================================

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


//------------------------------------------------------------
//%%%%------ GET TRACKS AND LOOK AT DEDX OF START OF TRACK.
TVector3 BestTrackletPos;
TVector3 BestTrackletDir;
std::vector<int> UsedTracklets;
	//for(auto const& trk : * tracks){
	for(size_t track_index = 0; track_index < ev_track->size(); track_index++) {
	//Check to see if the track is not already accounted for
		bool usedtrk = false;
		for(unsigned int check=0; check<UsedTracklets.size();check++) if(track_index == UsedTracklets[check]) usedtrk=true;
		if(usedtrk) continue;
	//%%%%%--- Get the start of the track 
	std::cout<<"Trajectory Points "<< ev_track->at(track_index).NumberTrajectoryPoints()<<std::endl;
	std::cout<<"Total Track Length "<< ev_track->at(track_index).Length()<<std::endl;
	auto const& hit_index_v = track_hit_ass[track_index];
	std::cout<<" Size of hit vector "<<hit_index_v.size()<<std::endl;
	auto thehit = ev_hit->at( hit_index_v[0] ) ;
       //auto dedxval = fCaloAlg.dEdx_AMP(h, 3.0);// What does pitch mean.... isn;t it just 3 cm
	
	//DEDX avg 
	double dedxavg = 0;
	for(unsigned int b=0; b<hit_index_v.size();b++){
		auto thehit = ev_hit->at( hit_index_v[b] ) ;
		unsigned int theplane = thehit.View();
	       auto dedxval = fCaloAlg.dEdx_AMP(thehit.Integral(),(thehit.PeakTime()+ tick_offset)*geom->TimeToCm(), 3.0, theplane );// What does pitch mean.... isn;t it just 3 cm
	       dedxavg += fCaloAlg.dEdx_AMP(thehit.Integral(),(thehit.PeakTime()+ tick_offset)*geom->TimeToCm(), 3.0, theplane );// What does pitch mean.... isn;t it just 3 cm
//	std::cout<<"dEdx value :"<<dedxval<<std::endl;
		}
	dedxavg/=hit_index_v.size();
	std::cout<<"dEdx avg value :"<<dedxavg<<std::endl;

	//Log info for Best tracklet with dEdx
		for(size_t pt=0; pt< ev_track->at(track_index).NumberTrajectoryPoints(); pt++){
			 std::cout<<"\t Position at point ("<< pt << ")   : [" <<ev_track->at(track_index).LocationAtPoint(pt)[0]<<","<<ev_track->at(track_index).LocationAtPoint(pt)[1]<<" , " <<ev_track->at(track_index).LocationAtPoint(pt)[2]<<" ]"<<std::endl;
			 std::cout<<"\t Direction at point ("<< pt << ")   : [" <<ev_track->at(track_index).DirectionAtPoint(pt)[0]<<","<<ev_track->at(track_index).DirectionAtPoint(pt)[1]<<" , " <<ev_track->at(track_index).DirectionAtPoint(pt)[2]<<" ]"<<std::endl;
			/*
			//Shitty hack;;;;;;;
			 auto planes = larutil::Geometry::GetME()->Views();
			 std::vector<larlite::geo::View_t> cplanes;
			 std::copy(planes.begin(), planes.end(), std::back_inserter(cplanes));
			std::cout<<" Size of planes "<< cplanes.size()<<std::endl;
			 auto colplane = cplanes[2];
			std::cout<<"Break point A"<<std::endl;
			if(trk.NumberdQdx(colplane)!=0){
			//Shitty hack^^^^^^^;
			 std::cout<<"DQdx Plane 2 "<< trk.DQdxAtPoint(pt, colplane)<<std::endl;
			}
			*/
			 //std::cout<<"DEdx "<< DQdxAtPoint<<std::endl;
			}// Loop over all the points
		BestTrackletPos = ev_track->at(track_index).Vertex();
		BestTrackletDir = ev_track->at(track_index).VertexDirection();
		}// for loop over tracks
	
	
	

//%%%%------ Take the Best Track track and then use the direction to contain the hits inside a cone.
	std::cout<<"&&&&CONE STUFF&&& "<<coneintpc<<std::endl;
	auto coneintpc = fgeoconic.ConeInTPC(BestTrackletPos,BestTrackletDir,ConeLength,angle, smoothness);
	std::cout<<"Is the cone in the TPC? "<<coneintpc<<std::endl;
	//if(coneintpc){
	auto ConeEdge = fgeoconic.FitConicalFeatures(BestTrackletPos,BestTrackletDir,ConeLength,angle, 2,smoothness);
	std::cout<<std::endl;
		for(int ep=0; ep<ConeEdge.size(); ep++)
			std::cout<<ConeEdge[ep].w<<" "<<ConeEdge[ep].t+250.0<<std::endl;
	std::cout<<std::endl;
	//}
	
	//%%%% Check to see if it's profile is more shower like or track like... 
		//&&&& base this on profile of hit spread or base it on profile shape
	//%%%%------ Then we need to see about refining this axis.... if that is the case?  // do we even want to try this here? 
	
//%%%%------ Take the leftover tracks and look for dEdx that are not contained in the previous cone.
	//%%%%------ acounted charge/ tracks 


	
//------------------------------------------------------------











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
