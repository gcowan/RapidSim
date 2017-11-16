#ifndef RAPIDIPSMEAR_H
#define RAPIDIPSMEAR_H

class RapidIPSmear {
	public:
		virtual ~RapidIPSmear() {}

		virtual std::pair<double,double> smearIP(double ip,double pt)=0;
};

#endif
