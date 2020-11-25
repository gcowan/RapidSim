#ifndef RAPIDACCEPTANCELHCB_H
#define RAPIDACCEPTANCELHCB_H

#include <vector>

#include "TLorentzVector.h"

#include "RapidAcceptance.h"

class RapidAcceptanceLHCb : public RapidAcceptance {
	public:
		RapidAcceptanceLHCb(AcceptanceType type, const std::vector<RapidParticle*>& parts, const std::vector<RapidCut*>& cuts)
			: RapidAcceptance(type, parts, cuts),
			  zC_(5.265), zTracker_(8.6),
			  xSizeTracker_(8.6*0.35), ySizeTracker_(8.6*0.27),
			  xMinTracker_(7.5*0.012), yMinTracker_(7.5*0.012)
		{}

		virtual void getDefaultPtRange(double& min, double& max);
		virtual void getDefaultEtaRange(double& min, double& max);

	private:
		virtual bool partInAcceptance(RapidParticle* part);
		virtual bool partInDownstream(RapidParticle* part);

		TLorentzVector magnetKick(TLorentzVector& vec, double charge);

		//parameters for determining whether a track remains in acceptance after the magnet

		//distance from PV to the magnet centre
		const double zC_;

		//z-pos of the end of the tracker
		const double zTracker_;

		//transverse size of the tracker
		const double xSizeTracker_;
		const double ySizeTracker_;

		//inner edge of the tracker
		const double xMinTracker_;
		const double yMinTracker_;
};

#endif
