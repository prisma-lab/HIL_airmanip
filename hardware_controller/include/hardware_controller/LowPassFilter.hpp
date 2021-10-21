
#ifndef _LowPassFilter_hpp_
#define _LowPassFilter_hpp_

#include <cmath>

class LowPassFilter{
public:
	//constructors
	LowPassFilter();
	LowPassFilter(float iCutOffFrequency, float iDeltaTime);
	//functions
	float update(float input);
	float update(float input, float deltaTime, float cutoffFrequency);
	//get and configure funtions
	float getOutput() const{return output;}
	void reconfigureFilter(float deltaTime, float cutoffFrequency);
private:
	float output;
	float ePow;
};

class NotchFilter { //500Hz, 1.9Hz cutoff
	public:
		NotchFilter();
		double update(double u);
	private:
		double u1,u2,y1,y2;
};

#endif //_LowPassFilter_hpp_
