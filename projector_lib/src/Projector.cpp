#include <projector_lib/Projector.hpp>
#include <iostream>
#include <corecrt_math_defines.h>

// constructor

Projector::Projector(std::unique_ptr<float[]> inputImg, std::shared_ptr<Geometry> geometry, SumAlgorithm sumAlgorithm, 
	int imgSize_x, int imgSize_y, int imgSize_z)
	: imgSize_x(imgSize_x), imgSize_y(imgSize_y), imgSize_z(imgSize_z), inputImg(MyImg(std::move(inputImg), imgSize_x, imgSize_y, imgSize_z)), geometry(std::move(geometry)), sumAlgorithm(sumAlgorithm){};


// common methods

std::pair<Point, Point> Projector::getIntersectionPoints(const Line& line) const {
	// Method to get intersection points of line and image or pixel edges.
	// Default method will get intersection points of line and image edges.

	int count = 0;
	int min_x, min_y, max_x, max_y;

	min_x = 0;
	min_y = 0;
	if (line.transpose) {
		max_x = imgSize_y;
		max_y = imgSize_x;
	}
	else {
		max_x = imgSize_x;
		max_y = imgSize_y;
	}
	
	double x_bottom = line.coor(min_y);
	double y_left = line.value(min_x);
	double x_top= line.coor(max_y);
	double y_right = line.value(max_x);

	bool bottom = (x_bottom >= min_x && x_bottom <= max_x);
	bool left = (y_left >= min_y && y_left <= max_y);
	bool top = (x_top >= min_x && x_top <= max_x);
	bool right = (y_right >= min_y && y_right <= max_y);
	count += (int)top + (int)bottom + (int)right + (int)left;

	if (count > 2) {
		if (top && bottom) {
			left = false;
			right = false;
		}
		else if (left && right) {
			top = false;
			bottom = false;
		}
		count = 2;
	}
	
	if (count < 2)
		return std::pair<Point, Point>(Point(), Point());

	Point point1, point2;
	double i, j, i_max, j_max;
	if (bottom) {
		i = x_bottom;
		j = min_y;
	}
	if (left) {
		i = min_x;
		j = y_left;
	}
	if (top) {
		i_max = x_top;
		j_max = max_y;
	}
	if (right) {
		i_max = max_x;
		j_max = y_right;
	}

	point1 = Point(std::floor(i) - 1, std::floor(j) - 1);
	point2 = Point(std::ceil(i_max) + 1, std::ceil(j_max) + 1);
	return std::pair<Point, Point>(point1, point2);
}

Line Projector::constructLine(const Line& line) const {
	// constructor of the line to determine transformation needed to 
	// lead line to the form with k>1 and create such a line with a new parameters

	Line out_line = line;
	if (line.isHorizontal) {
		out_line.isVertical = true;
		out_line.isHorizontal = false;
		out_line.transpose = true;
		return out_line;
	}
	if (line.isVertical) { return out_line; }

	out_line.b = line.b + inputImg.center_x * line.k - inputImg.center_y;

	if (line.k < 0) {
		out_line.reverse_x = true;
		out_line.k = -out_line.k;
		out_line.angle = M_PI - out_line.angle;
	}
	if (std::abs(line.k) < 1) {
		out_line.transpose = true;
		out_line.b = out_line.coor(0.0);
		out_line.k = 1 / out_line.k;
		out_line.angle = M_PI_2 - out_line.angle;
	}
	if (!out_line.transpose) {
		out_line.b = out_line.b - inputImg.center_x * out_line.k + inputImg.center_y;
	}
	else {
		out_line.b = out_line.b - inputImg.center_y * out_line.k + inputImg.center_x;
	}
	out_line.reverse_k = 1.0 / out_line.k;
	out_line.reverse_ab = 1.0 / std::pow(1.0 + out_line.k * out_line.k, 0.5);
	return out_line;
}

Line Projector::constructLine(const Line& line, bool transpose, bool reverse_x) const {
	// overloaded constructor to create line with fixed transpose and reverse parameters

	Line out_line(line);
	out_line.transpose = transpose;
	out_line.reverse_x = reverse_x;
	if (line.isHorizontal && transpose) {
		out_line.isVertical = true;
		out_line.isHorizontal = false;
		return out_line;
	}
	if (line.isVertical && reverse_x) {
		out_line.b = imgSize_x - line.b; // 2a - x where a inverse line
		return out_line;
	}
	out_line.b = line.b + inputImg.center_x * line.k - inputImg.center_y;

	if (reverse_x) {
		out_line.k = -out_line.k;
		out_line.angle = M_PI - out_line.angle;
	}
	if (transpose) {
		out_line.b = out_line.coor(0.0);
		out_line.k = 1 / out_line.k;
		out_line.angle = M_PI_2 - out_line.angle;
	}
	if (!transpose) {
		out_line.b = out_line.b - inputImg.center_x * out_line.k + inputImg.center_y;
	}
	else {
		out_line.b = out_line.b - inputImg.center_y * out_line.k + inputImg.center_x;
	}

	out_line.reverse_k = 1.0 / out_line.k;
	out_line.reverse_ab = 1.0 / std::pow(1.0 + out_line.k * out_line.k, 0.5);

	return out_line;
}

Line Projector::constructLine3D(const Line& line) const {
	Line out_line = line;

	if (line.kx < 0) {
		out_line.reverse_x = true;
		out_line.kx = -line.kx;
		out_line.x0 = inputImg.size_x - line.x0;
	}
	if (line.ky < 0) {
		out_line.reverse_y = true;
		out_line.ky = -line.ky;
		out_line.y0 = inputImg.size_y - line.y0;
	}
	if (line.kz < 0) {
		out_line.reverse_z = true;
		out_line.kz = -line.kz;
		out_line.z0 = inputImg.size_z - line.z0;
	}
	double a, b, c;
	a = out_line.getLambda(0.0, 0);
	b = out_line.getLambda(0.0, 1);
	c = out_line.getLambda(0.0, 2);
	double lambda = std::max({ out_line.getLambda(0.0, 0), out_line.getLambda(0.0, 1), out_line.getLambda(0.0, 2) });
	Point p = out_line.getPoint(lambda);

	assert(p.x < INFINITY);
	assert(p.y < INFINITY);
	assert(p.z < INFINITY);

	out_line.x0 = p.x;
	out_line.y0 = p.y;
	out_line.z0 = p.z;

	return out_line;
}

std::pair<Line, Line> Projector::sortLines(const Line& line1, const Line& line2) const {
	// Return pair (upper, lower) lines. For vertical lines return (right, left) lines.

	Line line_hv, line_o;
	if (line1.isHorizontal || line1.isVertical) {
		line_hv = line1;
		line_o = line2;
	}
	else {
		line_hv = line2;
		line_o = line1;
	}

	bool transpose = false;
	bool reverse_x = false;
	if (line_hv.isHorizontal) {
		transpose = true;
	}
	if (line_o.k < 0) {
		reverse_x = true;
	}

	return std::pair<Line, Line>(constructLine(line_hv, transpose, reverse_x), constructLine(line_o, transpose, reverse_x));

}



double Projector::singlePixelArea(int pixel_i, int pixel_j, const Line& line) const {
	// return the part of the pixel above the line
	// (1 - sum) would denote the part below the line due to whole area being 1.0
	double sum = 0.0;
	double dist = (pixel_j + 0.5 - line.k * (pixel_i + 0.5) - line.b) * line.reverse_ab;
	if (std::abs(dist) > 0.8) {
		if (dist >= 0) {
			sum = 1.0;
		}
		else {
			sum = 0.0;
		}
		//return sum;
	}
	double x1 = line.coor(pixel_j);
	double x2 = line.coor((pixel_j + 1));
	double y1 = line.value(pixel_i);
	double y2 = line.value((pixel_i + 1));
	int count = 0;
	bool bottom = (x1 >= pixel_i && x1 < (pixel_i + 1));
	bool top	= (x2 > pixel_i && x2 <= (pixel_i + 1));
	bool left	= (y1 >= pixel_j && y1 < (pixel_j + 1));
	bool right	= (y2 > pixel_j && y2 <= (pixel_j + 1));
	count += (int)top + (int)bottom + (int)right + (int)left;
	if (count > 2) {
		if (top && bottom) {
			left = false;
			right = false;
		}
		else if (left && right) {
			top = false;
			bottom = false;
		}
		count = 2;
	}

	x1 -= std::floor(x1);
	x2 = 1.0 - (std::ceil(x2) - x2);
	y1 -= std::floor(y1);
	y2 = 1.0 - (std::ceil(y2) - y2);

	if (count < 2) {
		if (dist >= 0) {
			sum = 1.0;
		}
		else {
			sum = 0.0;
		}
		return sum;
	}
	if (line.k > 0) {
		if (bottom && top) {
			sum = (x2 + x1) / 2.0;
		}
		else if (left && right) {
			sum = 1.0 - (y2 + y1) / 2.0;
		}
		else if (bottom) {
			assert(right);
			sum = 1.0 - (1.0 - x1) * (y2) / 2.0;
		}
		else if (top) {
			assert(left);
			sum = (1.0 - y1) * (x2) / 2.0;
		}
	}
	else if (line.k < 0) {
		if (bottom && top) {
			sum = 1.0 - (x2 + x1) / 2.0;
		}
		else if (left && right) {
			sum = 1.0 - (y2 + y1) / 2.0;
		}
		else if (bottom) {
			assert(left);
			sum = 1.0 - (x1 * y1) / 2.0;
		}
		else if (top) {
			assert(right);
			sum = (1.0 - x2) * (1.0 - y2) / 2.0;
		}
	}
	
	assert(0.0 <= sum && sum <= 1.0);
	
	return sum;
}

float Projector::manyPixelArea(int i_min, int i_max, int j, bool upper, const Line& line) const {

	float sum = 0.0f;
	int size_x = line.transpose ? imgSize_y : imgSize_x;
	int size_y = line.transpose ? imgSize_x : imgSize_y;
	if (j < 0 || j >= size_y) {
		return sum;
	}
	if (i_min < 0)
		i_min = -1;
	if (i_max >= size_x)
		i_max = size_x;
	for (int i = i_min; i < i_max; i++) {
		float a = inputImg.get(i, j, line.transpose, line.reverse_x);
		if (a == 0.0f) {
			continue;
		}
		double b = singlePixelArea(i, j, line);
		if (!upper) {
			b = 1.0 - b;
		}
		sum += a * b;
	}
	return sum;	
}

// legacy
float Projector::sumAreaExact(const Line& line_1, const Line& line_2) const {
	// Summation algorithm for the area between two arbitrary lines.
	// One of the line always have k>1 the other can be vertical or with k<-1.

	float sum = 0.0f;
	float sum1, sum2, sum3;
	double x_left_1, x_right_1;
	double x_left_2, x_right_2;

	assert(!line_1.isHorizontal);
	if (line_1.isVertical && line_2.isVertical)
	{
		int s = line_1.transpose ? imgSize_x : imgSize_y;
		sum1 = 0;// sumNeibs(-1, s, std::floor(line_1.b), line_1.transpose, line_1.reverse_x);
		sum2 = 0;// sumNeibs(-1, s, std::floor(line_2.b), line_1.transpose, line_1.reverse_x);
		if (line_1.b > line_2.b) {
			return sum1 * (line_1.b - std::floor(line_1.b)) + sum2 * (1 - (line_2.b - std::floor(line_2.b)));
		}
		else {
			return sum1 * (1 - (line_1.b - std::floor(line_1.b))) + sum2 * (line_2.b - std::floor(line_2.b));
		}
	}

	std::pair<Point, Point> p_1, p_2;
	bool is_left = false;
	p_2 = getIntersectionPoints(line_2);
	if (line_1.isVertical) {
		p_1 = std::pair<Point, Point>(Point(line_1.b, -1), Point(line_1.b, imgSize_y + 1));
		is_left = line_1.b < line_2.coor(0.0);
	}
	else {
		p_1 = getIntersectionPoints(line_1);
	}

	bool sign_1 = line_1.k > 0, sign_2 = line_2.k > 0;

	int j, j_max;
	j = std::min({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });
	j_max = std::max({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });

	if (p_2.first.x == -999) {
		j = std::min(p_1.first.y, p_1.second.y);
	}
	assert(j_max >= j);
	if (line_1.isVertical) {
		x_left_1 = line_1.b;
		x_right_1 = line_1.b;
	}
	double left, right;
	
	while (true) {
		if (!line_1.isVertical) {
			x_left_1 = line_1.coor(j);
			x_right_1 = line_1.coor(j + 1);
		}
		x_left_2 = line_2.coor(j);
		x_right_2 = line_2.coor(j + 1);

		left = std::min({ x_left_1, x_right_1, x_left_2, x_right_2 });
		right = std::max({ x_left_1, x_right_1, x_left_2, x_right_2 });
		if (line_1.isVertical) {
			if (is_left) {
				sum2 = (line_1.b - std::floor(line_1.b)) * inputImg.get(std::floor(left), j, line_1.transpose, line_1.reverse_x);
				sum3 = manyPixelArea(std::floor(left), std::floor(right) + 1, j, sign_2, line_2);
			}
			else {
				sum2 = (1 - (line_1.b - std::floor(line_1.b))) * inputImg.get(std::floor(right), j, line_1.transpose, line_1.reverse_x);
				sum3 = manyPixelArea(std::floor(left), std::floor(right) + 1, j, !sign_2, line_2);
			}
			sum += std::abs(sum2 - sum3);
		}
		else {
			sum2 = manyPixelArea(std::floor(left), std::floor(right) + 1, j, sign_1, line_1);
			sum3 = manyPixelArea(std::floor(left), std::floor(right) + 1, j, sign_2, line_2);
			sum += std::abs(sum2 - sum3);
		}

		j += 1;
		if (j > j_max)
			return sum;
	}

}

// weight algos

void Projector::weightNeibsAreaDual(int i_min, int i_max, int j, const Line& line1, const Line& line2,
	int* coorDst, float* weightsDst, int* sizeDst) const {

	float sum = 0.0f;
	int size_x = line1.transpose ? imgSize_y : imgSize_x;
	int size_y = line1.transpose ? imgSize_x : imgSize_y;
	bool upper1 = line1.k > 0;
	bool upper2 = line2.k > 0;
	if (j < 0 || j >= size_y) {
		return;
	}
	if (i_min < 0)
		i_min = 0;
	if (i_max >= size_x)
		i_max = size_x;
	for (int i = i_min; i < i_max; i++) {
		int coor = inputImg.get_coor(i, j, line1.transpose, line1.reverse_x);
		if (coor < 0) {
			continue;
		}
		double area1 = singlePixelArea(i, j, line1);
		double area2 = singlePixelArea(i, j, line2);
		if (!upper1) {
			area1 = 1.0 - area1;
		}
		if (!upper2) {
			area2 = 1.0 - area2;
		}
		if (area1 == area2)
			continue;

		coorDst[*sizeDst] = coor;
		weightsDst[*sizeDst] = std::abs(area1-area2);
		*sizeDst += 1;
		
	}
	return;
}

void Projector::weightNeibsAreaSingle(int i_min, int i_max, int j, const Line& line, bool upper,
	int* coorDst, float* weightsDst, int* sizeDst) const {

	float sum = 0.0f;
	int size_x = line.transpose ? imgSize_y : imgSize_x;
	int size_y = line.transpose ? imgSize_x : imgSize_y;
	if (j < 0 || j >= size_y) {
		return;
	}
	if (i_min < 0)
		i_min = 0;
	if (i_max >= size_x)
		i_max = size_x;
	for (int i = i_min; i < i_max; i++) {
		int coor = inputImg.get_coor(i, j, line.transpose, line.reverse_x);
		if (coor < 0) {
			continue;
		}
		double area = singlePixelArea(i, j, line);
		if (!upper) {
			area = 1.0 - area;
		}
		coorDst[*sizeDst] = coor;
		weightsDst[*sizeDst] = area;
		*sizeDst += 1;

	}
	return;
}

void Projector::weightNeibsLine(double j_min, double j_max, double i, int* coorDst, float* weightsDst, int* sizeDst, bool transpose, bool reverse_x, int slice, float koeff) const {
	
	if (j_min < 0)
		j_min = 0;
	if (transpose) {
		if (j_max >= inputImg.size_x)
		{
			j_max = inputImg.size_x;
		}
	}
	else {
		if (j_max >= inputImg.size_y)
		{
			j_max = inputImg.size_y;
		}
	}
	for (int j = j_min; j < j_max; j++) {
		int coor = inputImg.get_coor(i, j, transpose, reverse_x, false, slice);
		if (coor >= 0) {
			coorDst[*sizeDst] = coor;
			weightsDst[*sizeDst] = koeff;
			*sizeDst += 1;
		}
	}

}

void Projector::getWeightsLine(const Line& line, int* coorDst, float* weightsDst, int* sizeDst, int slice, bool isBinary)  const {
	// main method for summation along the line based on the distance that the beam
	// travels inside the pixel

	if (line.isVertical) {
		int s = line.transpose ? imgSize_x : imgSize_y;
		weightNeibsLine(-1, s, line.b, coorDst, weightsDst, sizeDst, line.transpose, line.reverse_x, slice);
		return;
	}

	std::pair<Point, Point> intersect = getIntersectionPoints(line);

	if (intersect.first.x == -999)
		return;

	int i, j, i_max, j_max;
	i = intersect.first.x;
	j = intersect.first.y;
	i_max = intersect.second.x;
	j_max = intersect.second.y;

	double new_x, new_y, weight;
	int new_i, new_j;
	assert(line.k >= 1.0);
	float temp_sum;
	int coor;
	while (true) {
		new_i = i + 1;
		new_y = line.value(new_i);
		new_j = std::floor(new_y);
		weightNeibsLine(j, new_j, i, coorDst, weightsDst, sizeDst, line.transpose, line.reverse_x, slice);
		if (isBinary) {
			weight = 0.5;
		}
		else {
			weight = std::abs(new_y - new_j);
		}
		coor = inputImg.get_coor(i, new_j, line.transpose, line.reverse_x, false, slice);
		if (coor >= 0) {
			coorDst[*sizeDst] = coor;
			weightsDst[*sizeDst] = weight;
			*sizeDst += 1;
		}
		coor = inputImg.get_coor(new_i, new_j, line.transpose, line.reverse_x, false, slice);
		if (coor >= 0) {
			coorDst[*sizeDst] = coor;
			weightsDst[*sizeDst] = (1.0f - weight);
			*sizeDst += 1;
		}

		i = new_i;
		j = new_j + 1;

		if (j > j_max) {
			break;
		}
	}
	float koeff = 1 / std::abs(std::sin(line.angle));
	for (int i = 0; i < *sizeDst; ++i) {
		weightsDst[i] *= koeff;
	}
	return;
}

void Projector::getWeightsArea(const Line& line_1, const Line& line_2, int* coorDst, float* weightsDst, int* sizeDst) const {

	float sum = 0.0f;
	float sum1, sum2, sum3;
	double x_left_1, x_right_1;
	double x_left_2, x_right_2;

	assert(!line_1.isHorizontal);
	if (line_1.isVertical && line_2.isVertical)
	{
		int s = line_1.transpose ? imgSize_x : imgSize_y;
		float koeff_1, koeff_2;
		if (line_1.b > line_2.b) {
			koeff_1 = (line_1.b - std::floor(line_1.b)); 
			koeff_2 = (1 - (line_2.b - std::floor(line_2.b)));
		}
		else {
			koeff_1 = (1 - (line_1.b - std::floor(line_1.b)));
			koeff_2 = (line_2.b - std::floor(line_2.b));
		}
		weightNeibsLine(-1, s, std::floor(line_1.b), coorDst, weightsDst, sizeDst, line_1.transpose, line_1.reverse_x, 0, koeff_1);
		weightNeibsLine(-1, s, std::floor(line_2.b), coorDst, weightsDst, sizeDst, line_1.transpose, line_1.reverse_x, 0, koeff_2);
		return;		
	}

	std::pair<Point, Point> p_1, p_2;
	p_1 = getIntersectionPoints(line_1);
	p_2 = getIntersectionPoints(line_2);

	bool sign_1 = line_1.k > 0, sign_2 = line_2.k > 0;

	int j, j_max;
	j = std::min({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });
	j_max = std::max({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });

	if (p_2.first.x == -999) {
		j = std::min(p_1.first.y, p_1.second.y);
	}
	assert(j_max >= j);
	double left, right;
	while (true) {
		x_left_1 = line_1.coor(j);
		x_right_1 = line_1.coor(j + 1);
		x_left_2 = line_2.coor(j);
		x_right_2 = line_2.coor(j + 1);

		left = std::min({ x_left_1, x_right_1, x_left_2, x_right_2 });
		right = std::max({ x_left_1, x_right_1, x_left_2, x_right_2 });

		weightNeibsAreaDual(std::floor(left), std::floor(right) + 1, j, line_1, line_2, coorDst, weightsDst, sizeDst);

		j += 1;
		if (j > j_max)
			return;
	}
}

void Projector::getWeightsAreaExact(const Line& line_1, const Line& line_2, int* coorDst, float* weightsDst, int* sizeDst) const {

	double x_left_1, x_right_1;
	double x_left_2, x_right_2;

	assert(!line_1.isHorizontal);
	if (line_1.isVertical && line_2.isVertical)
	{
		int s = line_1.transpose ? imgSize_x : imgSize_y;
		float koeff_1, koeff_2;
		if (line_1.b > line_2.b) {
			koeff_1 = (line_1.b - std::floor(line_1.b));
			koeff_2 = (1 - (line_2.b - std::floor(line_2.b)));
		}
		else {
			koeff_1 = (1 - (line_1.b - std::floor(line_1.b)));
			koeff_2 = (line_2.b - std::floor(line_2.b));
		}
		weightNeibsLine(-1, s, std::floor(line_1.b), coorDst, weightsDst, sizeDst, line_1.transpose, line_1.reverse_x, 0, koeff_1);
		weightNeibsLine(-1, s, std::floor(line_2.b), coorDst, weightsDst, sizeDst, line_1.transpose, line_1.reverse_x, 0, koeff_2);
		return;
	}


	std::pair<Point, Point> p_1, p_2;
	bool is_left = false;
	p_2 = getIntersectionPoints(line_2);
	if (line_1.isVertical) {
		p_1 = std::pair<Point, Point>(Point(line_1.b, -1), Point(line_1.b, imgSize_y + 1));
		is_left = line_1.b < line_2.coor(0.0);
	}
	else {
		p_1 = getIntersectionPoints(line_1);
	}

	bool sign_1 = line_1.k > 0, sign_2 = line_2.k > 0;

	int j, j_max;
	j = std::min({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });
	j_max = std::max({ p_1.first.y, p_1.second.y, p_2.first.y, p_2.second.y });

	if (p_2.first.x == -999) {
		j = std::min(p_1.first.y, p_1.second.y);
	}
	if (line_1.isVertical) {
		x_left_1 = line_1.b;
		x_right_1 = line_1.b;
	}
	assert(j_max >= j);
	double left, right;

	
	bool upper = !is_left ^ sign_2;
	while (true) {
		if (!line_1.isVertical) {
			x_left_1 = line_1.coor(j);
			x_right_1 = line_1.coor(j + 1);
		}
		x_left_2 = line_2.coor(j);
		x_right_2 = line_2.coor(j + 1);

		left = std::min({ x_left_1, x_right_1, x_left_2, x_right_2 });
		right = std::max({ x_left_1, x_right_1, x_left_2, x_right_2 });
		if (line_1.isVertical) {
			int coor;
			float value;

			if (is_left) {
				coor = inputImg.get_coor(std::floor(left), j, line_1.transpose, line_1.reverse_x);
				value = - (line_1.b - std::floor(line_1.b));
			}
			else {
				coor = inputImg.get_coor(std::floor(right), j, line_1.transpose, line_1.reverse_x);
				value = - (1 - (line_1.b - std::floor(line_1.b)));
			}
			coorDst[*sizeDst] = coor;
			weightsDst[*sizeDst] = value;
			*sizeDst += 1;
			weightNeibsAreaSingle(std::floor(left), std::floor(right) + 1, j, line_2, upper, coorDst, weightsDst, sizeDst);
		}
		else {
			weightNeibsAreaDual(std::floor(left), std::floor(right) + 1, j, line_1, line_2, coorDst, weightsDst, sizeDst);
		}
		j += 1;
		if (j > j_max)
			return;
	}
}

void Projector::getWeightsLine3D(const Line& line, int* coorDst, float* weightsDst, int* sizeDst, bool isBinary) const {

	double lambda, lambda_min = 0, lambda_max = std::min({ line.getLambda(inputImg.size_x, 0), line.getLambda(inputImg.size_y, 1), line.getLambda(inputImg.size_z, 2) });
	lambda = lambda_min;

	Point p = line.getPoint(lambda);
	std::vector<double> lambdas(3);
	int coors[3] = {};
	coors[0] = std::floor(p.x), coors[1] = std::floor(p.y), coors[2] = std::floor(p.z);

	for (int i = 0; i < 3; ++i) {
		lambdas[i] = line.getLambda(coors[i] + 1.0, i);
	}
	int index = std::min_element(lambdas.begin(), lambdas.end()) - lambdas.begin();
	float weight = 1;
	while (true) {
		if (!isBinary) {
			weight = (lambdas[index] - lambda);
		}
		int target_coor = inputImg.get_coor(coors[0], coors[1], coors[2], line.reverse_x, line.reverse_y, line.reverse_z);
		if (target_coor >= 0) {
			coorDst[*sizeDst] = target_coor;
			weightsDst[*sizeDst] = weight;
			*sizeDst += 1;
		}
		
		lambda = lambdas[index];
		coors[index] += 1;
		lambdas[index] = line.getLambda(coors[index] + 1.0, index);
		index = std::min_element(lambdas.begin(), lambdas.end()) - lambdas.begin();
		assert(lambdas[index] >= lambda);
		
		if (lambda > lambda_max) {
			return;
		}
	}
}


// main methods

std::unique_ptr<unsigned char[]> Projector::getForwardProjectionImage() {
	// normalize and return full projection in an image format
	// x axis is angle to detector as an index in input 'angles' vector
	// y axis is a detector pixel

	if (forwardProjection == nullptr)
		buildForwardProjection();

	std::unique_ptr<unsigned char[]> outCharImg(new unsigned char[geometry->nDetectors * geometry->angles.size()]);
	float max_val = *std::max_element(&forwardProjection[0], &forwardProjection[geometry->nDetectors * geometry->angles.size() - 1]);
	// check for too small max value; if (max < 1.0) then max = 1.0
	float reverse_max = (max_val > 1.0f) ? (1 / max_val) : 1.0f;
	for (int i = 0; i < geometry->angles.size(); ++i) {
		for (int j = 0; j < geometry->nDetectors; ++j) {
			outCharImg[i * geometry->nDetectors + j] = (unsigned char)(forwardProjection[i * geometry->nDetectors + j] * 255.0f * reverse_max);
		}
	}
	return outCharImg;
};

std::unique_ptr<float[]> Projector::getForwardProjection() {
	// getter for python interface

	if (!forwardProjection) {
		buildForwardProjection();
	}
	return std::move(forwardProjection);
}

// getWeightsAreaExact(line1, line2, coorDst, weightsDst, sizeDst);
// main methods weights

void Projector::buildBackProjection() {
	// build back projection
	int sizeDst = 0;
	int weights_size = imgSize_x + imgSize_y + imgSize_z + 1;
	if (sumAlgorithm == SumAlgorithm::AREA || sumAlgorithm == SumAlgorithm::AREA_EXACT)
		weights_size *= 2;
	int* coorDst = new int[weights_size]();
	float* weightsDst = new float[weights_size]();

	int sizeZ = geometry->imgCenterZ > 0 ? geometry->imgCenterZ * 2 : 1;
	
	backProjection = std::unique_ptr<float[]>(new float[geometry->imgCenterX * geometry->imgCenterY * sizeZ * 4]{});
	for (int angleIndex = 0; angleIndex < geometry->angles.size(); ++angleIndex) {
		for (int detectorIndex = 0; detectorIndex < geometry->nDetectors; ++detectorIndex) {
			getWeights(angleIndex, detectorIndex, coorDst, weightsDst, &sizeDst);
			float value = 0.0f;
			if (sizeDst > 0) {
				for (int i = 0; i < sizeDst; ++i)
					value += weightsDst[i];

				value = forwardProjection[angleIndex * geometry->nDetectors + detectorIndex] / value;
				if (value > 1) {
					float a = 1;
				}
				projectBack(value, coorDst, weightsDst, sizeDst);
			}
			
		}
	}
	delete[] coorDst;
	delete[] weightsDst;
}

void Projector::buildForwardProjection() {
	int sizeDst = 0;
	int weights_size = imgSize_x + imgSize_y + imgSize_z + 1;
	if (sumAlgorithm == SumAlgorithm::AREA || sumAlgorithm == SumAlgorithm::AREA_EXACT)
		weights_size *= 2;
	int* coorDst = new int[weights_size]();
	float* weightsDst = new float[weights_size]();
	forwardProjection = std::unique_ptr<float[]>(new float[geometry->nDetectors * geometry->angles.size()]{});

	for (int angleIndex = 0; angleIndex < geometry->angles.size(); ++angleIndex) {
		for (int detectorIndex = 0; detectorIndex < geometry->nDetectors; ++detectorIndex) {
			
			getWeights(angleIndex, detectorIndex, coorDst, weightsDst, &sizeDst);
			projectForward(angleIndex, detectorIndex, coorDst, weightsDst, sizeDst);
		}
	}

	delete[] coorDst;
	delete[] weightsDst;
}

void Projector::projectBack(float value, int* coors, float* weights, int size) {
	for (int i = 0; i < size; ++i) {
		backProjection[coors[i]] += weights[i] * value;
	}
}

void Projector::projectForward(int angleIndex, int detectorIndex, int* coors, float* weights, int size) {
	float sum = 0.0f;
	for (int i = 0; i < size; ++i) {
		sum += weights[i] * inputImg.inputImg[coors[i]];
	}
	forwardProjection[angleIndex * geometry->nDetectors + detectorIndex] = sum;
}

void Projector::getWeights(int angleIndex, int detectorIndex, int* coorDst, float* weightsDst, int* sizeDst) {

	Line line, line1, line2;
	std::pair<Line, Line> line_pair;
	int detector_i, detector_j;
	int slice;
	*sizeDst = 0;
	switch (sumAlgorithm)
	{
	case SumAlgorithm::LINE:
		line = constructLine(geometry->getLine(angleIndex, detectorIndex));
		getWeightsLine(line, coorDst, weightsDst, sizeDst);
		break;
	case SumAlgorithm::LINE3DPARALLEL:
		detector_i = detectorIndex % geometry->nDetectorsX;
		detector_j = detectorIndex / geometry->nDetectorsX;
		slice = std::floor((detector_j + 0.5) - geometry->nDetectorsY*0.5 + geometry->imgCenterZ);
		line = constructLine(geometry->getLine(angleIndex, detector_i));
		getWeightsLine(line, coorDst, weightsDst, sizeDst, slice);
		break;
	case SumAlgorithm::LINE3DFAN:
		line = constructLine3D(geometry->getLine(angleIndex, detectorIndex));
		getWeightsLine3D(line, coorDst, weightsDst, sizeDst);
		break;
	case SumAlgorithm::BINARY:
		line = constructLine(geometry->getLine(angleIndex, detectorIndex));
		getWeightsLine(line, coorDst, weightsDst, sizeDst, 0, true);
		break;
	case SumAlgorithm::BINARY3DPARALLEL:
		detector_i = detectorIndex % geometry->nDetectorsX;
		detector_j = detectorIndex / geometry->nDetectorsX;
		slice = std::floor((detector_j + 0.5) - geometry->nDetectorsY * 0.5 + geometry->imgCenterZ);
		line = constructLine(geometry->getLine(angleIndex, detector_i));
		getWeightsLine(line, coorDst, weightsDst, sizeDst, slice, true);
		break;
	case SumAlgorithm::BINARY3DCONE:
		line = constructLine3D(geometry->getLine(angleIndex, detectorIndex));
		getWeightsLine3D(line, coorDst, weightsDst, sizeDst, true);
		break;
	case SumAlgorithm::AREA:
		line_pair = geometry->getLinePair(angleIndex, detectorIndex);
		line1 = constructLine(line_pair.first);
		line2 = constructLine(line_pair.second, line1.transpose, line1.reverse_x);
		getWeightsArea(line1, line2, coorDst, weightsDst, sizeDst);
		break;
	case SumAlgorithm::AREA_EXACT:
		line_pair = geometry->getLinePair(angleIndex, detectorIndex, false);
		if ((line_pair.first.isHorizontal ^ line_pair.second.isHorizontal) || (line_pair.first.isVertical ^ line_pair.second.isVertical)) {
			std::pair<Line, Line> sorted_lines = sortLines(line_pair.first, line_pair.second);
			line1 = sorted_lines.first;
			line2 = sorted_lines.second;
		}
		else {
			line1 = constructLine(line_pair.first);
			line2 = constructLine(line_pair.second, line1.transpose, line1.reverse_x);
		}
		getWeightsAreaExact(line1, line2, coorDst, weightsDst, sizeDst);
		break;

	default:
		break;
	}
}
