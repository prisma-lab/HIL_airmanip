
#include "../include/hardware_controller/LowPassFilter.hpp"

#define ERROR_CHECK (true)

#if ERROR_CHECK
#include <iostream>
#endif

LowPassFilter::LowPassFilter():
	output(0),
	ePow(0){}

LowPassFilter::LowPassFilter(float iCutOffFrequency, float iDeltaTime):
	output(0),
	ePow(1-exp(-iDeltaTime * 2 * M_PI * iCutOffFrequency))
{
	#if ERROR_CHECK
	if (iDeltaTime <= 0){
		std::cout << "Warning: A LowPassFilter instance has been configured with 0 s as delta time.";
		ePow = 0;
	}
	if(iCutOffFrequency <= 0){
		std::cout << "Warning: A LowPassFilter instance has been configured with 0 Hz as cut-off frequency.";
		ePow = 0;
	}
	#endif
}

float LowPassFilter::update(float input){
	return output += (input - output) * ePow;
}

float LowPassFilter::update(float input, float deltaTime, float cutoffFrequency){
	reconfigureFilter(deltaTime, cutoffFrequency); //Changes ePow accordingly.
	return output += (input - output) * ePow;
}

void LowPassFilter::reconfigureFilter(float deltaTime, float cutoffFrequency){
	#if ERROR_CHECK
	if (deltaTime <= 0){
		std::cout << "Warning: A LowPassFilter instance has been configured with 0 s as delta time.";
		ePow = 0;
	}
	if(cutoffFrequency <= 0){
		std::cout << "Warning: A LowPassFilter instance has been configured with 0 Hz as cut-off frequency.";
		ePow = 0;
	}
	#endif
	ePow = 1-exp(-deltaTime * 2 * M_PI * cutoffFrequency);
}


NotchFilter::NotchFilter() {
	u1=u2=y1=y2=0;
}


double NotchFilter::update(double u) {
	double Num[] = {0.979482760981449, -1.958810850280416, 0.979482760981449};
	double Den[] = {-1.958810850280416, 0.958965521962899};
	//Numerator   = [0.9882 -1.9758 0.9882];  % Numerator coefficient vector
    //Denominator = [1 -1.9758 0.9764];        % Denominator coefficient vector

    //vecu = [u u1 u2];
    //vecy = [y1 y2];
    //y = (Numerator*vecu')-(Denominator(2:end)*vecy');
	double y = Num[0]*u + Num[1]*u1 + Num[2]*u2 - Den[0]*y1 - Den[1]*y2;
    
    y2 = y1;
    y1 = y;
    u2 = u1;
    u1 = u;

	return y;
}