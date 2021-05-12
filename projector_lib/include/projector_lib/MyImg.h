#ifndef _MYIMG_H_
#define _MYIMG_H_
#include <memory>

class MyImg {
public:
	const int size_x, size_y, size_z;
	const double center_x, center_y, center_z;
	std::unique_ptr<float[]> inputImg;

	MyImg(std::unique_ptr<float[]> inputImg, int size_x, int size_y, int size_z = 1);


	float get(int i, int j, bool transpose = false, bool reverse_x = false, bool reverse_y = false, int slice = 0) const;
	int get_coor(int i, int j, bool transpose = false, bool reverse_x = false, bool reverse_y = false, int slice = 0) const;
	float get(int i, int j, int k, bool reverse_x = false, bool reverse_y = false, bool reverse_z = false) const;
	int get_coor(int i, int j, int k, bool reverse_x, bool reverse_y, bool reverse_z) const;
};

#endif