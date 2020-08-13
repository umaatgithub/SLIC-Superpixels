#include <opencv2/highgui.hpp>
#include "slic6D_implementation/slic6d.h"

#include <iostream>
using namespace std;

int main(int argc, char **argv){
//    if (argc < 3) {
//        cout << "Usage: metrics_eval <input_dir> <outputFile>";
//        return 1;
//    }


    string image_path = "/home/ura/Documents/KUL/Projects/Vansteelandt/Data/Ortho/ORTHO stad_0_0.png";
    string dsm_path = "/home/ura/Documents/KUL/Projects/Vansteelandt/Data/DSM/DSM_stad_25831_NN_0_0.tif";
    //int nr_superpixels = atoi(argv[2]);

    cv::Mat image_ortho = cv::imread(image_path);
    cv::imshow( "Input image", image_ortho );
    cv::waitKey();
    std::cout << "Image type : " << image_ortho.type() << std::endl;

    cv::Vec3b pixel = image_ortho.at<cv::Vec3b>(0,0);
    std::cout << "Image pixel : " << (int)pixel.val[0] << " " << (int)pixel.val[1] << " " << (int)pixel.val[2] << std::endl;

    cv::Mat image_dsm = cv::imread(dsm_path, cv::IMREAD_UNCHANGED);
    cv::imshow("Input DSM image", image_dsm);
    cv::waitKey();
    double minVal, maxVal;
    cv::minMaxLoc(image_dsm, &minVal, &maxVal);
    std::cout << "DSM image min-max values : " << minVal << "     " << maxVal << std::endl;

    Slic6D slic;
    slic.generate_superpixels(&image_ortho, &image_dsm, 100);
    slic.display_clusters("Initial cluster", false);
    slic.create_connectivity();
    slic.display_clusters("Connected cluster", true);
}
