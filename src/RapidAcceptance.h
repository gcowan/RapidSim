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

		enum DetectorType {
			LHCB
		};

		static RapidAcceptance::AcceptanceType typeFromString(TString str);

		static RapidAcceptance::DetectorType detectorFromString(TString str);

		RapidAcceptance(AcceptanceType type, const std::vector<RapidParticle*>& parts, const std::vector<RapidCut*>& cuts)
			: type_(type),
			  cuts_(cuts)
		{setup(parts);}

        virtual ~RapidAcceptance() {}

		virtual bool isSelected();

	private:
		void setup(std::vector<RapidParticle*> parts);

		virtual bool inAcceptance();

		virtual bool motherInAcceptance();
		virtual bool allInAcceptance();
		virtual bool allInDownstream();

		virtual bool partInAcceptance(RapidParticle* part)=0;
		virtual bool partInDownstream(RapidParticle* part)=0;

		AcceptanceType type_;

		RapidParticle* top_;

		std::vector<RapidParticle*> parts_;

		std::vector<RapidCut*> cuts_;
};

#endif
