#include "CImg.h"
#include <projector_lib/Projector.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <projector_lib/Geometries.h>
#include <chrono>

std::vector<double> getAngles(double start, double end, int count) {
    std::vector<double> angles;
    double delta = (end - start) / count;
    for (int i = 0; i < count; ++i){
        angles.push_back(start + delta*i);
    }
    return angles;
}

std::string getPath(std::string type, std::string dim, std::string file_type, bool out = false, std::string model = "", std::string algo = "", std::string geometry = "") {
    std::string out_string = "D:/Projects/tomo_cmake/tomo_cmake/projector_test/";
    if (out) {
        out_string += "out_images/";
        if (type.size() > 0)
            out_string += type;
        if (dim.size() > 0)
            out_string += "_" + dim;
        out_string += "." + file_type;
    }
    else {
        out_string += "assets/";
        if (type == "sinogram") {
            out_string += "target_sinograms/";
        }
        if (type == "phantom") {
            out_string += "phantoms/"+model;
        }
        else if (type.size() > 0)
            out_string += type;
        if (dim.size() > 0)
            out_string += "_" + dim;
        if (type != "phantom") {
            if (model.size() > 0)
                out_string += "_" + model;
        }        
        if (geometry.size() > 0)
            out_string += "_" + geometry;
        if (algo.size() > 0) {
            if (algo == "binary") {
                out_string += "_line";
            }
            else {
                out_string += "_" + algo;
            }
        }

        out_string += "." + file_type;
    }
    return out_string;
}

std::unique_ptr<Projector> getProjector(std::string dim, int angle_count = 180, int detector_count_x = 256, int detector_count_y = 1
    , std::string model = "", std::string algo_string = "", std::string geometry_string = "", double distance_source = 7000.0, double distance_obj = 500.0, std::unique_ptr<float[]> input_img = nullptr
    , int size_x = 1, int size_y =1, int size_z = 1) {
    bool exact = false;
    if (algo_string == "fan")
        bool exact = true;
    std::string type = "bmp";
    if (dim == "3D") {
        type = "txt";
    }
    bool custom_phantom = false;
    std::string p = getPath("phantom", dim, type, false, model);
    std::vector<double> angles = getAngles(0, M_PI, angle_count);
    std::shared_ptr<Geometry> input_geometry = nullptr;
    auto sum_algo = SumAlgorithm::LINE;
    std::unique_ptr<float[]> img;
    if (input_img) {
        img = std::move(input_img);
        custom_phantom = true;
    }
    if (dim == "2D"){
        if (!custom_phantom) {
            const char* const phantom_path = p.c_str();
            cimg_library::CImg<unsigned char> phantom(phantom_path);
            size_x = phantom.width();
            size_y = phantom.height();
            img = std::unique_ptr<float[]>(new float[size_x * size_y]{});
            for (int i = 0; i < size_x; ++i) {
                for (int j = 0; j < size_y; ++j) {
                    img[i + j * size_x] = (float)phantom(i, j, 0, 0);
                }
            }

        }
        if (algo_string == "area") {
            if (exact)
                sum_algo = SumAlgorithm::AREA_EXACT;
            else
                sum_algo = SumAlgorithm::AREA;
        }
        else if (algo_string == "line") {
            sum_algo = SumAlgorithm::LINE;
        }
        if (geometry_string == "fan") {
            input_geometry = std::move(std::make_shared<GeometryFanBeam>(GeometryFanBeam(angles, detector_count_x, 1.0, distance_source, distance_obj, size_x * 0.5, size_y * 0.5)));
        }
        else {
            input_geometry = std::move(std::make_shared<GeometryParallel>(GeometryParallel(angles, detector_count_x, 1.0, size_x * 0.5, size_y * 0.5)));
        }

    }

    if (dim == "3D") {
        if (!custom_phantom) {
            std::ifstream input_phantom(p);
            // read phantom
            input_phantom >> size_x >> size_y >> size_z;
            img = std::unique_ptr<float[]>(new float[size_x * size_y * size_z]{});
            for (int k = 0; k < size_z; ++k) {
                for (int j = 0; j < size_y; ++j) {
                    for (int i = 0; i < size_x; ++i) {
                        float a;
                        input_phantom >> a;
                        img[i + size_x * j + size_x * size_y * k] = a;
                    }
                }
            }
            input_phantom.close();
        }
        if (geometry_string == "fan") {
            input_geometry = std::move(std::make_shared<GeometryFanBeam3D>(GeometryFanBeam3D(angles, detector_count_x, detector_count_y, 1.0, distance_source, distance_obj, size_x * 0.5, size_y * 0.5, size_z * 0.5)));
            sum_algo = SumAlgorithm::LINE3DFAN;
        }
        if (geometry_string == "par") {
            input_geometry = std::move(std::make_shared<GeometryParallel3D>(GeometryParallel3D(angles, detector_count_x, detector_count_y, 1.0, size_x * 0.5, size_y * 0.5, size_z * 0.5)));
            sum_algo = SumAlgorithm::LINE3DPARALLEL;
        }
    }

    return std::unique_ptr<Projector>(new Projector(std::move(img), std::move(input_geometry), sum_algo, size_x, size_y, size_z));
}


void customLineProjections() {
    std::vector<double> angles = getAngles(0.0, 180.0, 5.0);
    int angle_count = angles.size();
    int detector_size = 256;
    int size_x = 256, size_y = 256;
    for (int k = 1; k <= 14; ++k) {
        std::string s = std::to_string(k);
        std::ifstream input_f("phantom\\"+s+".txt");
        std::unique_ptr<float[]> img(new float[size_x * size_y]);

        unsigned char color[] = { 255 };
        cimg_library::CImg<float> cimage(size_x, size_y, 1, 1, 0);

        for (int i = 0; i < size_x; ++i) {
            for (int j = 0; j < size_y; ++j) {
                float a;
                input_f >> a;
                img[i * size_x + j] = a;
                //color[0] = a*255;
                //cimage.draw_point(i, j, color);
            }
        }

        input_f.close();
        Projector myProj(std::move(img), std::make_shared<GeometryParallel>(GeometryParallel(angles, detector_size, 1.0, size_x * 0.5, size_y * 0.5))
            , SumAlgorithm::LINE, size_x, size_y);
        std::ofstream output_f("phantom\\out_sino\\" + s + "_small.txt");
        //std::ofstream output_f("phantom\\out_sino\\" + s + "_large.txt");
        output_f.setf(std::ios::scientific);
        output_f.precision(std::numeric_limits<double>::digits10 + 1);
        for (int i = 0; i < angle_count; ++i) {
            for (int j = 0; j < detector_size; ++j) {
                Line line = myProj.geometry->getLine(angles[i], j);
                float sum = myProj.sumLine(line);
                float x1 = 0.0f, x2 = 256.0f;
                float y1 = line.value(x1), y2 = line.value(x2);
                if (std::abs(line.k) > 1) {
                    y1 = 0.0f;
                    y2 = 256.0f;
                    x1 = line.coor(y1);
                    x2 = line.coor(y2);
                }
                x1 -= size_x * 0.5f;
                x2 -= size_x * 0.5f;
                y1 -= size_y * 0.5f;
                y2 -= size_y * 0.5f;
                if (line.isVertical) {
                    x1 = line.b - size_x * 0.5f;
                    x2 = line.b - size_x * 0.5f;
                    y1 =        - size_y * 0.5f;
                    y2 =          size_y * 0.5f;
                }
                if (line.isHorizontal) {
                    x1 =        - size_x * 0.5f;
                    x2 =          size_x * 0.5f;
                    y1 = line.b - size_y * 0.5f;
                    y2 = line.b - size_y * 0.5f;
                }
                output_f << x1 << " " << y1 << " " << x2 << " " << y2 << " " << sum << std::endl;
            }
        }
        output_f.close();
        //cimage.display();
        std::cerr << k << " finished" << std::endl;

    }
    return;
}

void test3D(std::string model) {
    double distance_obj = 500.0, distance_source = 7000.0;
    float delta_limit = 2000000.0;

    int size_x, size_y, size_z;
    int detector_size_x, detector_size_y, angle_count;
    std::string p = getPath("phantom", "3D", "txt", false, model);
    std::string s = getPath("sinogram", "3D", "txt", false, model);
    std::ifstream input_phantom(p);
    std::ifstream input_sinogram(s);
    std::string out_sino = getPath("out_sinogram", "", "bmp", true);
    std::string out_diff = getPath("out_residual_float", "", "bmp", true);
    const char* const out_sino_path = out_sino.c_str();
    const char* const out_diff_path = out_diff.c_str();
    
    // read phantom and target
    input_phantom >> size_x >> size_y >> size_z;
    std::unique_ptr<float[]> img(new float[size_x * size_y * size_z]{});
    for (int k = 0; k < size_z; ++k) {
        for (int j = 0; j < size_y; ++j) {
            for (int i = 0; i < size_x; ++i) {
                float a;
                input_phantom >> a;
                img[i + size_x * j + size_x * size_y * k] = a;
            }
        }
    }
    input_phantom.close();
    input_sinogram >> angle_count >> detector_size_x >> detector_size_y;
    std::unique_ptr<float[]> target(new float[detector_size_x * detector_size_y * angle_count]{});
    for (int k = 0; k < angle_count; ++k) {
        for (int j = 0; j < detector_size_y; ++j) {
            for (int i = 0; i < detector_size_x; ++i) {
                float a;
                input_sinogram >> a;
                target[i + detector_size_x * j + detector_size_y * detector_size_x * k] = a;
            }
        }
    }
    input_sinogram.close();
    std::vector<double> angles = getAngles(0.0, M_PI, angle_count);
    std::shared_ptr<Geometry> input_geometry = nullptr;
    input_geometry = std::move(std::make_shared<GeometryFanBeam3D>(GeometryFanBeam3D(angles, detector_size_x, detector_size_y, 1.0, distance_source, distance_obj, size_x * 0.5, size_y * 0.5, size_z * 0.5)));
    Projector proj(std::move(img), input_geometry, SumAlgorithm::LINE3DFAN, size_x, size_y, size_z);
    proj.buildForwardProjection();

    // count sigma and max delta
    float delta = 0.0f;
    float total_delta = 0.0f;
    float out_delta = 1.0f, max_delta = 0.0f;
    float max_val = 0.0f;
    int detector_count = detector_size_x * detector_size_y;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_count; ++j) {
            delta = target[i * detector_count + j] - proj.forwardProjection[i * detector_count + j];
            if (proj.forwardProjection[i * detector_count + j] > max_val) {
                max_val = proj.forwardProjection[i * detector_count + j];
            }
            if (std::abs(delta) > max_delta) {
                max_delta = std::abs(delta);
                out_delta = std::abs(delta) / proj.forwardProjection[i * detector_count + j];
            }

            float delta_sq = delta * delta;
            total_delta += delta_sq;
            if (delta == 0)
                continue;
                        
            if (std::abs(delta) > delta_limit) {
                std::cerr << "delta: " << delta << " target: " << target[i * detector_count + j] << " res: "
                    << proj.forwardProjection[i * detector_count + j] << " det: " << j << " ang: " << i << std::endl;

                proj.getLineProjectionTest(i, j);
            }
        }
    }
    float sigma = std::pow(total_delta / (angle_count * detector_count - 1), 0.5);
    std::cout << "sigma: " << sigma << " max_delta: " << max_delta << " delta/value: " << out_delta << std::endl;

    // draw sino and residual picture
    unsigned char color3[] = { 255 , 0, 0 };
    unsigned char color[] = {0.0f};
    cimg_library::CImg<unsigned char> image(detector_size_x, angle_count*detector_size_y, 1, 1, 0);
    cimg_library::CImg<unsigned char> diff_float(detector_size_x, angle_count * detector_size_y, 1, 3, 0);

    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_count; ++j) {
            int j_x = j % detector_size_x;
            int j_y = j / detector_size_x;

            float d = (target[i * detector_count + j] - proj.forwardProjection[i * detector_count + j]);
            float value = d / max_delta;
            int k = value > 0 ? 0 : 1;
            color3[0] = 0;
            color3[1] = 0;
            color3[k] = abs(value * 255);
            diff_float.draw_point(j_x, j_y + i * detector_size_y, color3);
            color[0] = proj.forwardProjection[i * detector_count + j]/max_val * 255;
            image.draw_point(j_x, j_y + i * detector_size_y, color);
        }
    }
    diff_float.normalize(0, 255);
    diff_float.save(out_diff_path);
    image.normalize(0, 255);
    image.save(out_sino_path);

}

void fullTest2D() {

    double delta_limit = 1000000;
    double delta_limit_2 = 2000000;
    double distance_obj = 1000.0, distance_source = 7000.0;

    std::string model = "phantom"; // "phantom", "model";
    std::string algo_string = "line";    // "area", "line", "binary";
    std::string geometry_string = "par";     // "fan", "par";
    int size_z = 1;
    int detector_size_y = 100;
    bool exact = true;
    bool by_weights = true;

    //std::string p = "..\\..\\..\\phantom" + model + geometry_string + algo_string + ".bmp";// ??
    std::string p = getPath("phantom", "2D", "bmp", false, model);
    std::string t = getPath("sinogram", "2D", "bmp", false, model, algo_string, geometry_string);
    std::string target_txt_path = getPath("sinogram", "2D", "txt", false, model, algo_string, geometry_string);

    std::string out_sino = getPath("out_sinogram", "", "bmp", true);
    std::string out_sino_target = getPath("out_sinogram_target", "", "bmp", true);
    std::string out_diff_float = getPath("out_residual_float", "", "bmp", true);
    std::string out_diff = getPath("out_residual", "", "bmp", true);
    std::string out_delta_txt_path = getPath("delta", "", "txt", true);
    std::string out_delta_symetry_txt_path = getPath("symetry_delta", "", "txt", true);

    const char* const phantom_path = p.c_str();
    const char* const target_path = t.c_str();
    const char* const out_sino_target_path = out_sino_target.c_str();
    const char* const out_sino_path = out_sino.c_str();
    const char* const out_diff_float_path = out_diff_float.c_str();
    const char* const out_diff_path = out_diff.c_str();

    cimg_library::CImg<unsigned char> phantom(phantom_path), target(target_path);
    std::ifstream input_f(target_txt_path);

    int size_x = phantom.width(), size_y = phantom.height();
    int detector_size = target.height(), angle_count = target.width();
    std::unique_ptr<float[]> img(new float[size_x * size_y]{});
    std::cerr << model << " " << geometry_string << " " << algo_string << " " << size_x << " " << size_y << " " << size_z << std::endl;

    // copy first channel of image to array assuming grayscale bmp
    for (int i = 0; i < size_x; ++i) {
        for (int j = 0; j < size_y; ++j) {
            img[i + j * size_x] = (float)phantom(i, j, 0, 0);
        }
    }
    //std::ifstream input_angles("angles.txt");
    //std::vector<double> angles;
    //input_angles >> angle_count;
    //for (int i = 0; i < angle_count; ++i) {
    //    double a;
    //    input_angles >> a;
    //    angles.push_back(a);
    //}

    std::vector<double> angles;
    for (int i = 0; i < angle_count; ++i) {
        angles.push_back((double)(i)*M_PI / 180.0);
    }

    auto sum_algo = SumAlgorithm::LINE;
    if (algo_string == "area") {
        if (exact)
            sum_algo = SumAlgorithm::AREA_EXACT;
        else
            sum_algo = SumAlgorithm::AREA;
    }
    else if (algo_string == "line") {
        if (by_weights) {
            sum_algo = SumAlgorithm::LINE_BY_WEIGHTS;
        }
        else {
            sum_algo = SumAlgorithm::LINE;
        }
    }
    else if (algo_string == "binary") {
        sum_algo = SumAlgorithm::BINARY;
    }
    std::shared_ptr<Geometry> input_geometry = nullptr;
    if (geometry_string == "fan") {
        input_geometry = std::move(std::make_shared<GeometryFanBeam>(GeometryFanBeam(angles, detector_size, 1.0, distance_source, distance_obj, size_x * 0.5, size_y * 0.5)));
    }
    else {
        input_geometry = std::move(std::make_shared<GeometryParallel>(GeometryParallel(angles, detector_size, 1.0, size_x * 0.5, size_y * 0.5)));
    }

    Projector myProj(std::move(img), std::move(input_geometry), sum_algo, size_x, size_y);

    std::unique_ptr<unsigned char[]> fp = myProj.getForwardProjectionImage();

    // compare the differences in image form,
    // count non-zero differences in pixels for statistics
    float delta = 0.0f, summ = 0.0f;
    int count = 0;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            int a = (int)target(i, j, 0, 0);
            int b = (int)fp[i * detector_size + (j)];
            delta += std::abs(b - a);
            if (a != b) {
                count += 1;
            }
        }
    }
    if (count == 0)
        std::cout << "Perfect!" << std::endl;
    else
        std::cout << "total delta: " << delta << " count: " << count << " average: " << delta / count << std::endl;

    // build residual image
    unsigned char color[] = { 255 };
    unsigned char color2[] = { 255 , 0, 0 };
    cimg_library::CImg<unsigned char> image(angle_count, detector_size, 1, 1, 0);
    cimg_library::CImg<unsigned char> diff(angle_count, detector_size, 1, 3, 0);

    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            int b = fp[i * detector_size + j];
            color[0] = b;
            image.draw_point(i, j, color);
            int k = b > target(i, j, 0, 0) ? 0 : 1;
            color2[0] = 0;
            color2[1] = 0;
            color2[k] = abs(b - target(i, j, 0, 0));
            diff.draw_point(i, j, color2);
        }
    }

    // save normalized
    image.normalize(0, 255);
    image.save(out_sino_path);
    diff.normalize(0, 255);
    diff.save(out_diff_path);

    // read ground truth sinogram from file

    std::unique_ptr<float[]> target_f(new float[angle_count * detector_size]);
    float max_f = 0;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            float a;
            input_f >> a;
            target_f[i * detector_size + j] = a;
            if (a > max_f) {
                max_f = a;
            }
        }
    }

    cimg_library::CImg<float> target_sino_img(angle_count, detector_size, 1, 1, 0);
    float color_f[] = { 255.0f };
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            float a = target_f[i * detector_size + j];
            color_f[0] = a / max_f;
            target_sino_img.draw_point(i, j, color_f);

        }
    }
    target_sino_img.normalize(0, 255);
    target_sino_img.save(out_sino_target_path);
    input_f.close();

    // calculate sigma from ground truth
    delta = 0.0f;
    float total_delta = 0.0f;
    float out_delta = 1.0f, max_delta = 0.0f;
    int pos_delta = 0, neg_delta = 0;
    float percentage_delta = 0.0f, max_percentage_delta = 0.0f;
    std::ofstream output_f(out_delta_txt_path);
    output_f << "delta " << "target " << "res" << std::endl;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            delta = target_f[i * detector_size + j] - myProj.forwardProjection[i * detector_size + j];
            if (target_f[i * detector_size + j] > 0.1 && myProj.forwardProjection[i * detector_size + j] > 0.1) {
                percentage_delta = delta / std::max(target_f[i * detector_size + j], myProj.forwardProjection[i * detector_size + j]);
                if (std::abs(percentage_delta) > std::abs(max_percentage_delta)) {
                    max_percentage_delta = percentage_delta;

                }
            }
            if (std::abs(delta) > max_delta) {
                max_delta = std::abs(delta);
                out_delta = std::abs(delta) / myProj.forwardProjection[i * detector_size + j];
            }
            float delta_sq = delta * delta;
            total_delta += delta_sq;
            if (delta == 0)
                continue;
            if (delta > 0)
                pos_delta += 1;
            else
                neg_delta += 1;
            if (std::abs(delta) > delta_limit_2 && i < angle_count / 2) {
                float delta1 = target_f[i * detector_size + j] - target_f[(180 + i) * detector_size + j];
                float delta2 = myProj.forwardProjection[i * detector_size + j] - myProj.forwardProjection[(180 + i) * detector_size + j];
                std::cout << "delta: " << delta << " target: " << delta1 << " res: " << delta2 << std::endl;
            }
            if (std::abs(delta) > 0) {
                Line line = myProj.geometry->getLine(myProj.geometry->angles[i], j);

                output_f << delta << " " << target_f[i * detector_size + j] << " "
                    << myProj.forwardProjection[i * detector_size + j] << std::endl;
            }
            if (std::abs(delta) > delta_limit) {
                std::cerr << "delta: " << delta << " delta_p: " << percentage_delta << " target: " << target_f[i * detector_size + j] << " res: "
                    << myProj.forwardProjection[i * detector_size + j] << " det: " << j << " ang: " << i << std::endl;

                myProj.getLineProjectionTest(i, j);
            }
        }
    }
    float sigma = std::pow(total_delta / (angle_count * detector_size - 1), 0.5);
    std::cout << "sigma: " << sigma << " max_delta: " << max_delta << " delta/value: " << out_delta << std::endl;
    std::cout << "pos: " << pos_delta << " neg: " << neg_delta << std::endl;
    output_f.close();

    // symmetry
    std::ofstream output_f1(out_delta_symetry_txt_path);
    output_f1 << "target " << "target_sym " << "res " << "res_sym " << "k " << "line_angle " << "angle_i " << "detector " << "delta " << "delta_sym_t " << "delta_sym_r " << std::endl;
    float t_delta_sym, r_delta_sym;
    for (int i = 0; i < angle_count / 2; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            Line line;
            float t_direct = target_f[i * detector_size + j];
            float t_sym = target_f[(180 + i) * detector_size + detector_size - j - 1];
            float r_direct = myProj.forwardProjection[i * detector_size + j];
            float r_sym = myProj.forwardProjection[(i + 180) * detector_size + detector_size - j - 1];
            delta = std::max(std::abs(t_direct - r_direct), std::abs(t_sym - r_sym));
            t_delta_sym = std::abs(t_direct - t_sym);
            r_delta_sym = std::abs(r_direct - r_sym);
            if (t_delta_sym < 1e-10 && r_delta_sym < 1e-10) {
                //continue;
            }
            int detector, angle_i;
            if (std::abs(t_direct - r_direct) > std::abs(t_sym - r_sym)) {
                line = myProj.geometry->getLine(i, j);
                detector = j;
                angle_i = i;
            }
            else {
                line = myProj.geometry->getLine((i + 180), detector_size - j - 1);
                detector = detector_size - j - 1;
                angle_i = 180 + i;
            }

            if (delta > 0) {
                output_f1 << t_direct << " " << t_sym << " " << r_direct << " " << r_sym << " " << line.k << " "
                    << line.angle << " " << angle_i << " " << detector << " " << delta << " " << t_delta_sym << " " << r_delta_sym << std::endl;
            }
        }
    }
    output_f1.close();

    // calculate residual picture from float values
    unsigned char color3[] = { 255 , 0, 0 };
    cimg_library::CImg<unsigned char> diff_float(angle_count, detector_size, 1, 3, 0);
    cimg_library::CImg<unsigned char> diff_percentage(angle_count, detector_size, 1, 3, 0);
    float max_p = 0.1f;
    for (int i = 0; i < angle_count; ++i) {
        for (int j = 0; j < detector_size; ++j) {
            float d = (target_f[i * detector_size + j] - myProj.forwardProjection[i * detector_size + j]);
            float value = d / max_delta;
            int k = value > 0 ? 0 : 1;
            color3[0] = 0;
            color3[1] = 0;
            color3[k] = abs(value * 255);
            diff_float.draw_point(i, j, color3);
            percentage_delta = d / std::max(target_f[i * detector_size + j], myProj.forwardProjection[i * detector_size + j]);
            percentage_delta = std::min(std::abs(percentage_delta), max_p);
            color3[k] = abs(percentage_delta / max_p * 255);
            diff_percentage.draw_point(i, j, color3);
        }
    }
    diff_float.normalize(0, 255);
    diff_float.save(out_diff_float_path);

    //diff_percentage.normalize(0, 255);
    //diff_percentage.save("out_residual_percentage.bmp");
}

void test3D_self_projection(){
    float delta_limit = 100000;

    int detector_count_x = 256;
    int detector_count_y = 21;
    int angle_count = 90;
    double distance_source  = 7000.0;
    double distance_obj     = 500.0;
    std::string model           = "model";
    std::string algo_string     = "line";
    std::string geometry_string = "par";
    
    std::string p = getPath("phantom", "2D", "txt", false, model);
    std::ifstream input_phantom(p);
    // read phantom
    int size_x = 1, size_y = 1, size_z = 1;
    input_phantom >> size_x >> size_y >> size_z;
    std::cerr << model << " " << geometry_string << " " << algo_string << " " << size_x << " " << size_y << " " << size_z << std::endl;
    std::unique_ptr<float[]> phantom2D(new float[size_x * size_y * size_z]{});
    for (int k = 0; k < size_z; ++k) {
        for (int j = 0; j < size_y; ++j) {
            for (int i = 0; i < size_x; ++i) {
                float a;
                input_phantom >> a;
                phantom2D[i + size_x * j + size_x * size_y * k] = a;
            }
        }
    }
    input_phantom.close();

    size_z = 121;
    std::unique_ptr<float[]> phantom3D(new float[size_x * size_y * size_z]{});
    for (int k = 0; k < size_z; ++k) {
        for (int j = 0; j < size_y; ++j) {
            for (int i = 0; i < size_x; ++i) {
                phantom3D[i + size_x * j + size_x * size_y * k] = phantom2D[i + size_x * j];
            }
        }
    }

    auto proj2D = getProjector("2D", angle_count, detector_count_x, 1,                model, algo_string, geometry_string, distance_source, distance_obj, std::move(phantom2D), size_x, size_y, 1);
    auto proj3D = getProjector("3D", angle_count, detector_count_x, detector_count_y, "", "", geometry_string, distance_source, distance_obj, std::move(phantom3D), size_x, size_y, size_z);

    std::string out_diff_float = getPath("out_residual_float", "", "bmp", true);
    const char* const out_diff_float_path = out_diff_float.c_str();
    std::string out_sino = getPath("out_sinogram", "", "bmp", true);
    std::string out_sino_2D = getPath("out_sinogram_target", "", "bmp", true);

    const char* const out_sino_path = out_sino.c_str();
    const char* const out_sino_2D_path = out_sino_2D.c_str();

    proj2D->buildForwardProjection();
    proj3D->buildForwardProjection();

    float max_val = 0.0f;
    float total_max_delta = 0.0f;
    double koeff = 1.0;
    for (int k = 0; k < detector_count_y; ++k) {
        if (geometry_string == "fan") {
            double angle = std::atan((detector_count_y * 0.5 - k - 0.5) / (distance_source + distance_obj));
            koeff = std::abs(std::cos(angle));
        }
        float delta = 0.0f;
        float total_delta = 0.0f;
        float out_delta = 1.0f, max_delta = 0.0f;
        
        for (int i = 0; i < angle_count; ++i) {
            for (int j = 0; j < detector_count_x; ++j) {
                float target = proj2D->forwardProjection[j + i * detector_count_x];
                float res = (proj3D->forwardProjection[j + k * detector_count_x + i * detector_count_x * detector_count_y]) * koeff;
                delta = target - res;
                
                if (std::abs(delta) > max_delta) {
                    max_delta = std::abs(delta);
                    out_delta = std::abs(delta) / std::max(res, target);
                    total_max_delta = std::abs(delta);
                }
                if (res > max_val) {
                    max_val = res;
                }
                float delta_sq = delta * delta;
                total_delta += delta_sq;

                //if (delta == 0)
                //    continue;
                if (std::abs(delta) > delta_limit) {
                    std::cerr << "delta: " << delta << " target: " << target << " res: "
                        << res << " det: " << j << " ang: " << i << std::endl;
                    
                }
            }
        }
        float sigma = std::pow(total_delta / (angle_count * detector_count_x - 1), 0.5);
        std::cout << "row: " << k << " sigma: " << sigma << " max_delta: " << max_delta << " delta/value: " << out_delta << std::endl;
    }
    std::cerr << "max_val: " << max_val << " total_max_delta: " << total_max_delta << std::endl;
    unsigned char color[] = { 255 , 0, 0 };
    unsigned char color2[] = { 0 };
    cimg_library::CImg<unsigned char> diff(detector_count_y * angle_count, detector_count_x, 1, 3, 0);
    cimg_library::CImg<unsigned char> sino(detector_count_y * angle_count, detector_count_x, 1, 1, 0);
    cimg_library::CImg<unsigned char> sino_2D(detector_count_y * angle_count, detector_count_x, 1, 1, 0);
    for (int k = 0; k < detector_count_y; ++k) {
        if (geometry_string == "fan") {
            double angle = std::atan((detector_count_y * 0.5 - k - 0.5) / (distance_source + distance_obj));
            koeff = std::abs(std::cos(angle));
        }
        for (int i = 0; i < angle_count; ++i) {
            for (int j = 0; j < detector_count_x; ++j) {
                float target = proj2D->forwardProjection[j + i * detector_count_x];
                float res    = proj3D->forwardProjection[j + k * detector_count_x + i * detector_count_x * detector_count_y] * koeff;
                int g = res > target ? 0 : 1;
                color[0] = 0;
                color[1] = 0;
                color[g] = std::abs(res - target)/total_max_delta * 255;
                diff.draw_point(k + i * detector_count_y, j, color);
                color2[0] = res / max_val * 255;
                sino.draw_point(k + i * detector_count_y, j, color2);
                color2[0] = target / max_val * 255;
                sino_2D.draw_point(k + i * detector_count_y, j , color2);
            }
        }
    }

    diff.normalize(0, 255);
    diff.save(out_diff_float_path);
    sino.normalize(0, 255);
    sino.save(out_sino_path);
    sino_2D.normalize(0, 255);
    sino_2D.save(out_sino_2D_path);

}

void measure(std::string model, std::string algo_string, std::string geometry_string, std::string dim
    , int detector_count_x, int detector_count_y, int angle_count) {
    if (dim == "2D") {
        detector_count_y = 1;
    }
    std::unique_ptr<Projector> proj = getProjector(dim, angle_count, detector_count_x, detector_count_y, model, algo_string, geometry_string);
    auto start = std::chrono::high_resolution_clock::now();
    proj->buildForwardProjection();
    auto end = std::chrono::high_resolution_clock::now();
    auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cerr << dim << " " << geometry_string << " " << algo_string << std::endl;
    std::cerr << "img_size: " << proj->imgSize_x << " " << proj->imgSize_y << " " << proj->imgSize_z << std::endl;
    std::cerr << "detector_size: " << detector_count_x << " " << detector_count_y << std::endl;
    std::cerr << "angle_count: " << angle_count << std::endl;
    std::cerr << "time ms: " << time_ms << std::endl;
}   


int main()
{   
    //int detector_count_x = 128;
    //int detector_count_y = 128;
    //int angle_count = 180;
    //std::string model = "phantom";
    //std::string algo_string = "line";
    //std::string geometry_string = "fan";
    //std::string dim = "2D";

    //measure(model, algo_string, geometry_string, dim, detector_count_x, detector_count_y, angle_count);
    fullTest2D();
    //test3D("model");
    //test3D_self_projection();
}
