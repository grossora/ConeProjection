/**
 * \file ClusterShower.h
 *
 * \ingroup ClusterShower
 * 
 * \brief Class def header for a class ClusterShower
 *
 * @author ryan
 */

/** \addtogroup ClusterShower

    @{*/

#ifndef LARLITE_CLUSTERSHOWER_H
#define LARLITE_CLUSTERSHOWER_H

#include "Analysis/ana_base.h"
#include "LArUtil/LArUtilManager.h"
#include "LArUtil/PxUtils.h"
#include "LArUtil/LArUtilBase.h"
#include "LArUtil/LArUtil-TypeDef.h"
#include "LArUtil/TimeService.h"

#include "AnalysisAlg/CalorimetryAlg.h"

#include "DataFormat/DataFormat-TypeDef.h"
#include "DataFormat/hit.h"
#include "DataFormat/cluster.h"
#include "DataFormat/track.h"
#include "DataFormat/mcshower.h"
#include "DataFormat/data_base.h"

#include <math.h>       /* aTan */
#define PI 3.14159265358979323846264338327950288419
#include <stdlib.h>     /* srand, rand */

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"



#include "geoconic.h"
//#include "ctools.h"
//#include "conicalprojection.h"
//#include "coneprofile.h"


#include "ClusterRecoUtil/ClusterParamsAlg.h"
#include "ClusterRecoUtil/CRUHelper.h"
#include "EMShowerTools/EMShowerProfile.h"

namespace larlite {
  /**
     \class ClusterShower
     User custom analysis class made by SHELL_USER_NAME
   */
  class ClusterShower : public ana_base{
  
  public:

    /// Default constructor
    ClusterShower(){ _name="ClusterShower"; _fout=0;}

    /// Default destructor
    virtual ~ClusterShower(){}

    /** IMPLEMENT in ClusterShower.cc!
        Initialization method to be called before the analysis event loop.
    */ 
    virtual bool initialize();

    /** IMPLEMENT in ClusterShower.cc! 
        Analyze a data event-by-event  
    */
    virtual bool analyze(storage_manager* storage);

    /** IMPLEMENT in ClusterShower.cc! 
        Finalize method to be called after all events processed.
    */
    virtual bool finalize();

  protected:


        //::larlite::conicalprojection fconeproj;
        ::larlite::geoconic fgeoconic;
	//::calo::CalorimetryAlg *fCaloAlg;
	::calo::CalorimetryAlg fCaloAlg;

        //::larlite::ctools fctools;
        //::larlite::coneprofile fconeprofile;

//	double AxisLength = 200;// Setting this here... This is the length of the cone
        //double openingangle = 50.0; // magic number place holder for now. 
        double openingangle ; // magic number place holder for now. 
        int smoothness = 16;// would be nice if this was even... but this gives the smoothness of the edge of the polygon cone
	double energy= -999;
	//double angle = 60;// magic number place holder for now. 
	double angle = 40;// magic number place holder for now. 
                bool coneintpc = true;
                double ConeLength = 4*14;
	
// temp stuff that needs to be removed 
int Eventcounter = 0;

  };
}
#endif

//**************************************************************************
// 
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
