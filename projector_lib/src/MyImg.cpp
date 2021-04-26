#include <projector_lib/MyImg.h>

MyImg::MyImg(std::unique_ptr<float[]> i_inputImg, const int& i_size_x, const int& i_size_y)
	: size_x(i_size_x), size_y(i_size_y), inputImg(std::move(i_inputImg)), center_x(i_size_x * 0.5), center_y(i_size_y * 0.5)
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

float MyImg::get_c(const double& x, const double& y) const {
	if (x < -center_x || x >= center_x || y < -center_y || y >= center_y) {
		return 0;
	}
	return inputImg[floor(x + center_x) + floor(y + center_y) * size_x];
};