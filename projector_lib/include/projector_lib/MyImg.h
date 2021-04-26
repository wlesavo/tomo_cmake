#ifndef _MYIMG_H_
#define _MYIMG_H_
#include <memory>

class MyImg {
public:
	const int size_x, size_y;
	const double center_x, center_y;
	std::unique_ptr<float[]> inputImg;
	MyImg(std::unique_ptr<float[]> inputImg, const int& size_x, const int& size_y);
	float get(int i, int j, bool transpose = false, bool reverse_x = false, bool reverse_y = false) const;
	float get_c(const double& x, const double& y) const;
};

#endif