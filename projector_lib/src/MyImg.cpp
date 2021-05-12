#include <projector_lib/MyImg.h>
#include <iostream>

MyImg::MyImg(std::unique_ptr<float[]> inputImg, int size_x, int size_y, int size_z)
	: size_x(size_x), size_y(size_y), size_z(size_z), inputImg(std::move(inputImg))
	, center_x(size_x * 0.5), center_y(size_y * 0.5), center_z(size_z * 0.5)
{};

float MyImg::get(int i, int j, bool transpose, bool reverse_x, bool reverse_y, int slice) const {
	if (transpose) {
		int a = j;
		j = i;
		i = a;
	}
	if (i < 0 || i >= size_x || j < 0 || j >= size_y || slice < 0 || slice >= size_z) {
		return 0;
	}
	if (reverse_x)
		i = size_x - 1 - i;
	if (reverse_y)
		j = size_y - 1 - j;

	int coor;

	coor = i + j * size_x + slice * size_x * size_y;

	return inputImg[coor];
};


int MyImg::get_coor(int i, int j, bool transpose, bool reverse_x, bool reverse_y, int slice) const {
	if (transpose) {
		int a = j;
		j = i;
		i = a;
	}

	if (i < 0 || i >= size_x || j < 0 || j >= size_y || slice < 0 || slice >= size_z) {
		return -1;
	}
	if (reverse_x)
		i = size_x - 1 - i;
	if (reverse_y)
		j = size_y - 1 - j;

	int coor = i + j * size_x + slice * size_x * size_y;
	
	return coor;
}


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


int MyImg::get_coor(int i, int j, int k, bool reverse_x, bool reverse_y, bool reverse_z) const {
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
	return coor;
};
