#ifndef RAPIDACCEPTANCE_H
#define RAPIDACCEPTANCE_H

#include <vector>

#include "TLorentzVector.h"

class RapidCut;
class RapidParticle;

class RapidAcceptance {
	public:
		enum AcceptanceType {
			ANY,
			MOTHERIN,
			ALLIN,
			ALLDOWNSTREAM
		};

		static RapidAcceptance::AcceptanceType typeFromString(TString str);

		RapidAcceptance(AcceptanceType type, const std::vector<RapidParticle*>& parts, const std::vector<RapidCut*>& cuts)
			: type_(type),
			  cuts_(cuts),
			  zC_(5.4), ptkick_(1.2), zTracker_(9.5),
			  xSizeTracker_(9.5*0.3), ySizeTracker_(9.5*0.25),
			  xMinTracker_(9.5*0.001), yMinTracker_(9.5*0.001)
		{setup(parts);}

		bool isSelected();

	private:
		void setup(std::vector<RapidParticle*> parts);

		bool inAcceptance();

		bool motherInAcceptance();
		bool allInAcceptance();
		bool allInDownstream();

		bool partInAcceptance(RapidParticle* part);
		bool partInDownstream(RapidParticle* part);

		TLorentzVector magnetKick(TLorentzVector& vec, double charge);

		AcceptanceType type_;

		RapidParticle* top_;

		std::vector<RapidParticle*> parts_;

		std::vector<RapidCut*> cuts_;

		//parameters for determining whether a track remains in acceptance after the magnet

		//distance from PV to the magnet centre
		const double zC_;

		//kick from the magnet
		const double ptkick_;

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
