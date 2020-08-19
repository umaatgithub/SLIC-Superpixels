#ifndef SLIC6D_H
#define SLIC6D_H

#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>

#include <vector>
#include <iostream>

class Slic6D
{

private:
    cv::Mat *image;
    cv::Mat *dsm;
    std::vector<std::vector<float> > centers;
    cv::Mat *clusters;
    cv::Mat *conn_clusters;
    cv::Mat *distances;

    int step, num_superpixels;

    void init();
    /* Compute the distance between a center and an individual pixel. */
    float compute_dist(int center_index, cv::Vec3b color, cv::Point3f coordinate);
    /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
    cv::Point2i find_local_minimum(cv::Mat *image, cv::Point2i center);

    void setNum_superpixels(int value);
    void setImage(cv::Mat *value);
    void setDsm(cv::Mat *value);

public:
    Slic6D();
    ~Slic6D();

    void generate_superpixels(cv::Mat *input_image, cv::Mat *dsm, int num_superpixels);
    void create_connectivity();
    cv::Mat display_clusters(std::string win_name, bool connected);

    cv::Mat *getConn_clusters() const;
};

#endif // SLIC6D_H
