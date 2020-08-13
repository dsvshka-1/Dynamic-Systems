#pragma once

class DOPRI5
{
	DOPRI5() {};
	virtual ~DOPRI5() {};

public:
	template<typename T> 
	void sol();
};