#include <iostream>
#include <boost/filesystem.hpp>

#include <opencv2/highgui.hpp>
#include "slic6D_implementation/slic6d.h"

namespace fs = boost::filesystem;

int main(int argc, char **argv){

    fs::path image_dir("/home/ura/Documents/KUL/Projects/Vansteelandt/Data/Ortho/");
    fs::path dsm_dir("/home/ura/Documents/KUL/Projects/Vansteelandt/Data/DSM/");
    fs::path output_dir("/home/ura/Documents/KUL/Projects/Vansteelandt/Data/slic/");

    std::vector<fs::path> files;
    copy(fs::directory_iterator(image_dir), fs::directory_iterator(), back_inserter(files));
    sort(files.begin(), files.end());

    for (auto&& image_file : files) {

        std::string image_name = image_file.filename().string();
        std::cout << "Processing image : " << image_name << std::endl;
        cv::Mat image_ortho = cv::imread(image_dir.string()+image_name);
        if(image_ortho.cols < 1024 || image_ortho.rows < 1024){
            continue;
        }
        //cv::imshow( "Input image", image_ortho );
        //cv::waitKey();
        //std::cout << "Image type : " << image_ortho.type() << std::endl;

        //cv::Vec3b pixel = image_ortho.at<cv::Vec3b>(0,0);
        //std::cout << "Image pixel : " << (int)pixel.val[0] << " " << (int)pixel.val[1] << " " << (int)pixel.val[2] << std::endl;

        std::string dsm_name;
        dsm_name = "DSM_stad_25831_CUB" + image_name.substr(10,image_name.find_first_of(".")-10)+".tif";
        cv::Mat image_dsm = cv::imread(dsm_dir.string()+dsm_name, cv::IMREAD_UNCHANGED);
        if(image_dsm.cols < 1024 || image_dsm.rows < 1024){
            continue;
        }
        //cv::imshow("Input DSM image", image_dsm);
        //cv::waitKey();
        //double minVal, maxVal;
        //cv::minMaxLoc(image_dsm, &minVal, &maxVal);
        //std::cout << "DSM image min-max values : " << minVal << "     " << maxVal << std::endl;

        Slic6D slic;
        slic.generate_superpixels(&image_ortho, &image_dsm, 200);
        //slic.display_clusters("Initial cluster", false);
        slic.create_connectivity();
        cv::Mat output = slic.display_clusters("Connected cluster", true);
        std::string output_name;
        output_name = image_name.substr(0,image_name.find_first_of("."))+"_slic-cub.png";
        cv::imwrite(output_dir.string()+output_name, output);
    }

}
