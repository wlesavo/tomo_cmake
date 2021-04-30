#include <projector_lib/Projector.hpp>
#include <iostream>
#include <corecrt_math_defines.h>

// constructor

Projector::Projector(std::unique_ptr<float[]> i_inputImg, std::shared_ptr<Geometry> i_geometry, SumAlgo sumAlgorithm, 
	int i_imgSize_x, int i_imgSize_y, int i_imgSize_z)
	: imgSize_x(i_imgSize_x), imgSize_y(i_imgSize_y), imgSize_z(i_imgSize_z), inputImg(MyImg(std::move(i_inputImg), i_imgSize_x, i_imgSize_y, i_imgSize_z)), geometry(std::move(i_geometry)), sumAlgorithm(sumAlgorithm){};

// common methods

std::pair<Point, Point> Projector::getIntersectionPoints(const Line& line, int pixel_i, int pixel_j) const {
	// Method to get intersection points of line and image or pixel edges.
	// Default method will get intersection points of line and image edges.

	int count = 0;
	int min_x, min_y, max_x, max_y;
	// get intersection points

	if (pixel_i == -999) {
		min_x = 0;
		min_y = 0;
		max_x = imgSize_x;
		max_y = imgSize_y;
	}
	else {
		min_x = pixel_i;
		min_y = pixel_j;
		max_x = pixel_i + 1;
		max_y = pixel_j + 1;
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

	point1 = Point(std::floor(i)-1, std::floor(j)-1);
	point2 = Point(std::ceil(i_max)+1, std::ceil(j_max)+1);
	return std::pair<Point, Point>(point1, point2);
}

std::pair<Line, Line> Projector::sortLines(const Line& line1, const Line& line2) const {
	// return pair (upper, lower) lines. for vertical lines return (right, left) lines
	Line line_hv, line_o;
	if (line1.hor || line1.vert) {
		line_hv = line1; 
		line_o = line2;
	}
	else {
		line_hv = line2;
		line_o = line1;
	}

	bool transpose = false;
	bool reverse_x = false;
	if (line_hv.hor) {
		transpose = true;
	}
	if (line_o.k < 0) {
		reverse_x = true;
	}
	
	return std::pair<Line, Line>(constructLine(line_hv, transpose, reverse_x), constructLine(line_o, transpose, reverse_x));

}

float Projector::sumNeibs(double j_min, double j_max, double i, bool transpose, bool reverse_x) const {
	// summation along a vertical axis taking into account transpose and reverse transformations
	
	float summ = 0;

	if (j_min < 0)
		j_min = -1;
	if (transpose) {
		if (j_max >= inputImg.size_y)
		{
			j_max = inputImg.size_y;
		}
	}
	else {
		if (j_max >= inputImg.size_x)
		{
			j_max = inputImg.size_x;
		}
	}
	
	for (int j = j_min; j < j_max; j++) {
		summ += inputImg.get(i, j, transpose, reverse_x);
	}
	
	return summ;
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

	// todo: check for upper and lower pixel positions // fixed
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
	if (j < 0 || j >= inputImg.size_y) {
		return sum;
	}
	if (i_min < 0)
		i_min = -1;
	if (i_max >= inputImg.size_x)
		i_max = inputImg.size_x;
	for (int i = i_min; i < i_max; i++) {
		if (upper) {
			float a = inputImg.get(i, j, line.transpose, line.reverse_x);
			if (a == 0.0f) {
				continue;
			}
			double b = singlePixelArea(i, j, line);
			sum += a * b;
			//sum += singlePixelArea(i, j, line) * inputImg.get(i, j, line.transpose, line.reverse_x);
		}
		else if (!upper) {
			float a = inputImg.get(i, j, line.transpose, line.reverse_x);
			if (a == 0.0f) {
				continue;
			}
			double b = (1.0 - singlePixelArea(i, j, line));
			sum += a * b;
		}
	}
	return sum;	
}

Line Projector::constructLine(const Line& line) const {
	// constructor of the line to determine transformation needed to 
	// lead line to the form with k>1 and create such a line with a new parameters

	Line out_line = line;
	if (line.hor) {
		out_line.vert = true;
		out_line.hor = false;
		out_line.transpose = true;
		return out_line;
	}
	if (line.vert) { return out_line; }

	out_line.b = line.b + inputImg.center_x * line.k - inputImg.center_y;
	if (std::abs(line.k) < 1) {
		out_line.transpose = true;
		out_line.b = out_line.coor(0.0);
		out_line.k = 1 / out_line.k;
		out_line.angle = M_PI_2 - out_line.angle;
	}
	if (line.k < 0) {
		out_line.reverse_x = true;
		out_line.k = - out_line.k;
		out_line.angle = M_PI - out_line.angle;
	}
	out_line.b = out_line.b - inputImg.center_x * out_line.k + inputImg.center_y;

	out_line.reverse_k = 1.0 / out_line.k;
	out_line.reverse_ab = 1.0 / std::pow(1.0 + out_line.k * out_line.k, 0.5);
	return out_line;
}

Line Projector::constructLine(const Line& line, bool transpose, bool reverse_x) const {
	// overloaded constructor to create line with fixed transpose and reverse parameters

	Line out_line(line);
	out_line.transpose = transpose;
	out_line.reverse_x = reverse_x;
	if (line.hor && transpose) {
		out_line.vert = true;
		out_line.hor = false;
		return out_line;
	}
	if (line.vert && reverse_x) {
		out_line.b = imgSize_x - line.b; // 2a - x where a inverse line
		return out_line;
	}

	out_line.b = line.b + inputImg.center_x * line.k - inputImg.center_y;

	if (transpose) {
		out_line.b = out_line.coor(0.0);
		out_line.k = 1 / out_line.k;
		out_line.angle = M_PI_2 - out_line.angle;
	}
	if (reverse_x) {
		out_line.k = -out_line.k;
		out_line.angle = M_PI - out_line.angle;
	}

	out_line.b = out_line.b - inputImg.center_x * out_line.k + inputImg.center_y;

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
	double lambda = std::max({ out_line.getLambda(0.0, 0), out_line.getLambda(0.0, 1), out_line.getLambda(0.0, 2)});
	Point p = out_line.getPoint(lambda);

	assert(p.x < INFINITY);
	assert(p.y < INFINITY);
	assert(p.z < INFINITY);

	out_line.x0 = p.x;
	out_line.y0 = p.y;
	out_line.z0 = p.z;

	return out_line;
}
// sum algos

float Projector::sumLine(const Line& line) const {
	// main method for summation along the line based on the distance that the beam
	// travels inside the pixel


	if (line.vert)
		return sumNeibs(-1, imgSize_y, line.b, line.transpose, line.reverse_x);

	std::pair<Point, Point> intersect = getIntersectionPoints(line);

	if (intersect.first.x == -999)
		return 0.0f;

	float sum = 0.0f;
	// choose axis to sum along based on line parameters
	// axis = 1; dy = 0, dx = 1 || axis = 0; dy = 1, dx = 0
	
	int i, j, i_max, j_max;
	i = intersect.first.x;
	j = intersect.first.y;
	i_max = intersect.second.x;
	j_max = intersect.second.y;

	// main summation cycle
	// summation is performed by sequencually finding next intersections of line and 
	// vertical (axis = 0) or horizontal (axis = 1) pixel border, 
	// summing up to this point and adding proportional values for the intersection area
	double new_x, new_y, frac;
	int new_i, new_j;
	assert(line.k >= 1.0);
	while (true) {
		new_i = i + 1;
		new_y = line.value(new_i);
		new_j = std::floor(new_y);
		sum += sumNeibs(j, new_j, i, line.transpose, line.reverse_x);
		frac = std::abs(new_y - new_j);
		float b = inputImg.get(i, new_j, line.transpose, line.reverse_x) * frac + inputImg.get(new_i, new_j, line.transpose, line.reverse_x) * (1 - frac);
		if (b > 0) {
			float a = 0;
		}
		sum += b;
		i = new_i;
		j = new_j + 1;
		assert(sum >= 0);
		if (j > j_max)
			return sum / std::abs(std::sin(line.angle));
		
	}
}

float Projector::sumLinear(const Line& line) const {
	// not finished method for summation with counting pixel contributions
	// based on distance to neighbour pixels
	// DO NOT USE AS IS

	// todo
	if (line.hor)
		return 0.0f; //sumNeibs(-1, imgSize_x, line.b, 0);
	else if (line.vert)
		return 0.0f; //sumNeibs(-1, imgSize_y, line.b, 1);

	std::pair<Point, Point> intersect = getIntersectionPoints(line);

	if (intersect.first.x == -999)
		return 0.0f;

	float sum = 0.0f;
	// choose axis to sum along based on line parameters
	// axis = 1; dy = 0, dx = 1 || axis = 0; dy = 1, dx = 0
	int axis = (std::abs(line.k) >= 1.0);
	int sign = (line.k > 0) - (line.k < 0);
	int i, j, i_max, j_max;
	i = intersect.first.x;
	j = intersect.first.y;
	i_max = intersect.second.x;
	j_max = intersect.second.y;

	// main summation cycle
	// summation is performed by sequencually finding next intersections of line and 
	// vertical (axis = 0) or horizontal (axis = 1) pixel border, 
	// summing up to this point and adding proportional values for the intersection area
	double new_x, new_y, frac;
	int new_i, new_j;
	double left, right;
	while (true) {
		if (axis == 1) { // dy == 0
			new_x = line.coor(j + 0.5);
			left = std::floor(new_x - 0.5);
			right = left + 1;
			frac = std::abs(new_x - 0.5 - left);
			sum += inputImg.get(left, j) * (1 - frac) + inputImg.get(right, j) * frac;
			j = j + 1;
			assert(sum >= 0);
			if (j > j_max + 1)
				return sum / std::abs(std::sin(line.angle));
		}
		else if (axis == 0) { // dx == 0
			new_y = line.value(i + 0.5);
			left = std::floor(new_y - 0.5);
			right = left + 1;
			frac = std::abs(new_y - 0.5 - left);
			sum += inputImg.get(i, left) * (1 - frac) + inputImg.get(i, right) * frac;
			i = i + 1;
			assert(sum >= 0);
			if (i > i_max + 1)
				return sum / std::abs(std::cos(line.angle));
		}
	}
}

float Projector::sumArea(const Line& line_1, const Line& line_2) const {
	// summation algorithm for the area between two parallel lines
	// both lines always have k>1 or both vertical

	float sum = 0.0f;
	float sum1, sum2, sum3;
	double x_left_1, x_right_1;
	double x_left_2, x_right_2;

	assert(!line_1.hor);
	if (line_1.vert && line_2.vert)
	{
		sum1 = sumNeibs(-1, imgSize_y, std::floor(line_1.b), line_1.transpose, line_1.reverse_x);
		sum2 = sumNeibs(-1, imgSize_y, std::floor(line_2.b), line_1.transpose, line_1.reverse_x);
		if (line_1.b > line_2.b) {
			return sum1 * (line_1.b - std::floor(line_1.b)) + sum2 * (1 - (line_2.b - std::floor(line_2.b)));
		}
		else {
			return sum1 * (1 - (line_1.b - std::floor(line_1.b))) + sum2 * (line_2.b - std::floor(line_2.b));
		}
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
		
		sum2 = manyPixelArea(std::floor(left), std::floor(right) + 1, j, sign_1, line_1);
		sum3 = manyPixelArea(std::floor(left), std::floor(right) + 1, j, sign_2, line_2);
		sum += std::abs(sum2 - sum3);
		
		j += 1;
		if (j > j_max)
			return sum;
	}


}

float Projector::sumAreaExact(const Line& line_1, const Line& line_2) const {
	// summation algorithm for the area between two arbitrary lines
	// one of the line always have k>1 the other can be vertical or with k<-1

	float sum = 0.0f;
	float sum1, sum2, sum3;
	double x_left_1, x_right_1;
	double x_left_2, x_right_2;

	assert(!line_1.hor);
	if (line_1.vert && line_2.vert)
	{
		sum1 = sumNeibs(-1, imgSize_y, std::floor(line_1.b), line_1.transpose, line_1.reverse_x);
		sum2 = sumNeibs(-1, imgSize_y, std::floor(line_2.b), line_1.transpose, line_1.reverse_x);
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
	if (line_1.vert) {
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
	if (line_1.vert) {
		x_left_1 = line_1.b;
		x_right_1 = line_1.b;
	}
	double left, right;
	
	while (true) {
		if (!line_1.vert) {
			x_left_1 = line_1.coor(j);
			x_right_1 = line_1.coor(j + 1);
		}
		x_left_2 = line_2.coor(j);
		x_right_2 = line_2.coor(j + 1);

		left = std::min({ x_left_1, x_right_1, x_left_2, x_right_2 });
		right = std::max({ x_left_1, x_right_1, x_left_2, x_right_2 });
		if (line_1.vert) {
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

float Projector::sumLine3D(const Line& line) const {
	double lambda, lambda_min = 0, lambda_max = std::min({ line.getLambda(inputImg.size_x, 0), line.getLambda(inputImg.size_y, 1), line.getLambda(inputImg.size_z, 2) });
	lambda = lambda_min;
	Point p = line.getPoint(lambda);
	std::vector<double> lambdas(3);
	int coor[3] = {};
	coor[0] = std::floor(p.x), coor[1] = std::floor(p.y), coor[2] = std::floor(p.z);
	// todo: handle exact interseptions
	for (int i = 0; i < 3; ++i) {
		lambdas[i] = line.getLambda(coor[i] + 1.0, i);
	}
	int index = std::min_element(lambdas.begin(), lambdas.end()) - lambdas.begin();
	float sum = 0.0f;
	while (true) {
		float b = (lambdas[index] - lambda) * inputImg.get(coor[0], coor[1], coor[2], line.reverse_x, line.reverse_y, line.reverse_z);
		sum += b;
		lambda = lambdas[index];
		coor[index] += 1;
		lambdas[index] = line.getLambda(coor[index] + 1.0, index);
		index = std::min_element(lambdas.begin(), lambdas.end()) - lambdas.begin();
		assert(lambdas[index] >= lambda);
		
		if (lambda > lambda_max) {
			return sum;
		} 		
	}
	return sum;
}
// main methods

std::unique_ptr<float[]> Projector::getSingleProjection(int i_angle) const {
	// summation along particular angle of projection for every pixel of detector

	std::unique_ptr<float[]> outLineImg(new float[geometry->detectorCount]);
	Line line, line1, line2;
	std::pair<Line, Line> line_pair;
	
	for (int i = 0; i < geometry->detectorCount; ++i) {
		float out = 0;
		
		switch (sumAlgorithm)
		{
		case SumAlgo::LINE:
			line = constructLine(geometry->v_GetNextLineCenter(i_angle, i));
			out = sumLine(line);
			break;
		case SumAlgo::LINE3D:
			line = constructLine3D(geometry->v_GetNextLineCenter(i_angle, i));
			out = sumLine3D(line);
			break;
		case SumAlgo::LINEAR:
			line = constructLine(geometry->v_GetNextLineCenter(i_angle, i));
			out = sumLinear(line);
			break;
		case SumAlgo::BINARY:
			break;
		case SumAlgo::AREA:
			line_pair = geometry->v_GetNextLinePair(i_angle, i);
			line1 = constructLine(line_pair.first);
			line2 = constructLine(line_pair.second, line1.transpose, line1.reverse_x);
			out = sumArea(line1, line2);
			break;
		case SumAlgo::AREA_EXACT:
			line_pair = geometry->v_GetNextLinePair(i_angle, i, false);
			if ((line_pair.first.hor ^ line_pair.second.hor) || (line_pair.first.vert ^ line_pair.second.vert)) {
				std::pair<Line, Line> sorted_lines = sortLines(line_pair.first, line_pair.second);
				line1 = sorted_lines.first;
				line2 = sorted_lines.second;
			}
			else {
				line1 = constructLine(line_pair.first);
				line2 = constructLine(line_pair.second, line1.transpose, line1.reverse_x);
			}
			out = sumAreaExact(line1, line2);
			break;

		default:
			break;
		}

		outLineImg[i] = out;
	}

	return outLineImg;
}

void Projector::buildFullProjection(){
	// build projection for every angle
	// x axis is angle to detector as an index in input 'angles' vector
	// y axis is a detector pixel

	std::unique_ptr<float[]> outfloatImg(new float[geometry->detectorCount * geometry->angles.size()]);
	for (int i = 0; i < geometry->angles.size(); ++i) {
		std::unique_ptr<float[]> lineImg = getSingleProjection(i);
		for (int j = 0; j < geometry->detectorCount; ++j) {
			outfloatImg[i * geometry->detectorCount + j ] = lineImg[j];
		}
	}
	fullProjection = std::move(outfloatImg);
}

std::unique_ptr<unsigned char[]> Projector::getFullProjectionImage() {
	// normalize and return full projection in an image format
	// x axis is angle to detector as an index in input 'angles' vector
	// y axis is a detector pixel

	if (fullProjection == nullptr)
		buildFullProjection();

	std::unique_ptr<unsigned char[]> outCharImg(new unsigned char[geometry->detectorCount * geometry->angles.size()]);
	float max_val = *std::max_element(&fullProjection[0], &fullProjection[geometry->detectorCount * geometry->angles.size() - 1]);
	// check for too small max value; if (max < 1.0) then max = 1.0
	float reverse_max = (max_val > 1.0f) ? (1 / max_val) : 1.0f;
	for (int i = 0; i < geometry->angles.size(); ++i) {
		for (int j = 0; j < geometry->detectorCount; ++j) {
			outCharImg[i * geometry->detectorCount + j] = (unsigned char)(fullProjection[i * geometry->detectorCount + j] * 255.0f * reverse_max);
		}
	}
	return outCharImg;
};

std::unique_ptr<float[]> Projector::getFullProjection() {
	// getter for python interface

	if (!fullProjection) {
		buildFullProjection();
	}
	return std::move(fullProjection);
}

// debug methods

float Projector::getLineProjectionTest(const int& i_angle, const int& detector) {
	// test function for checking particular angle and detector

	float out = 0.0f;
	Line line, line1, line2;
	std::pair<Line, Line> line_pair;

	switch (sumAlgorithm)
	{
	case SumAlgo::LINE:
		line = constructLine(geometry->v_GetNextLineCenter(i_angle, detector));
		out = sumLine(line);
		break;
	case SumAlgo::LINE3D:
		line = constructLine3D(geometry->v_GetNextLineCenter(i_angle, detector));
		out = sumLine3D(line);
		break;
	case SumAlgo::LINEAR:
		line = constructLine(geometry->v_GetNextLineCenter(i_angle, detector));
		out = sumLinear(line);
		break;
	case SumAlgo::BINARY:
		break;
	case SumAlgo::AREA:
		line_pair = geometry->v_GetNextLinePair(i_angle, detector);
		line1 = constructLine(line_pair.first);
		line2 = constructLine(line_pair.second, line1.transpose, line1.reverse_x);
		out = sumArea(line1, line2);
		break;
	case SumAlgo::AREA_EXACT:
		line_pair = geometry->v_GetNextLinePair(i_angle, detector, false);
		if ((line_pair.first.hor ^ line_pair.second.hor) || (line_pair.first.vert ^ line_pair.second.vert)) {
			std::pair<Line, Line> sorted_lines = sortLines(line_pair.first, line_pair.second);
			line1 = sorted_lines.first;
			line2 = sorted_lines.second;
		}
		else {
			line1 = constructLine(line_pair.first);
			line2 = constructLine(line_pair.second, line1.transpose, line1.reverse_x);
		}
		out = sumAreaExact(line1, line2);
		break;
	default:
		break;
	}
	return out;
}
