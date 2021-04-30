#include <projector_lib/MyImg.h>

MyImg::MyImg(std::unique_ptr<float[]> i_inputImg, int i_size_x, int i_size_y, int i_size_z)
	: size_x(i_size_x), size_y(i_size_y), size_z(i_size_z), inputImg(std::move(i_inputImg))
	, center_x(i_size_x * 0.5), center_y(i_size_y * 0.5), center_z(i_size_z * 0.5)
{};

float MyImg::get(int i, int j, bool transpose, bool reverse_x, bool reverse_y) const{
	if (i < 0 || i >= size_x || j < 0 || j >= size_y) {
		return 0;
	}
	int coor;
	if (reverse_x)
		i = size_x - 1 - i;
	if (reverse_y)
		j = size_y - 1 - j;
	if (transpose)
		coor = j + i * size_y;
	else
		coor = i + j * size_x;
	return inputImg[coor];
};

float MyImg::get(int i, int j, int k, bool reverse_x, bool reverse_y, bool reverse_z) const {
	if (i < 0 || i >= size_x || j < 0 || j >= size_y || k < 0 || k >= size_z) {
		return 0;
	}
	int coor;
	if (reverse_x)
		i = size_x - 1 - i;
	if (reverse_y)
		j = size_y - 1 - j;
	if (reverse_z)
		k = size_z - 1 - k;
	coor = i + j * size_x + k * size_x * size_y;
	return inputImg[coor];
};
