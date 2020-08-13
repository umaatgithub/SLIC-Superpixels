#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>

using namespace std;

int main(){

    std::string path= "/home/ura/Documents/KUL/Projects/Vansteelandt/Data/DSM/DSM_stad_25831.tif";
    cv::Mat image_dsm = cv::imread(path, cv::IMREAD_UNCHANGED);
    cout << "Input image size : " << image_dsm.rows << " " << image_dsm.cols << endl;
    double minVal, maxVal;
    cv::minMaxLoc(image_dsm, &minVal, &maxVal);
    std::cout << "Input image min-max values : " << minVal << "     " << maxVal << std::endl;

    image_dsm.setTo(-1, image_dsm==-9999);
    cv::minMaxLoc(image_dsm, &minVal, &maxVal);
    std::cout << "Input min reset image min-max values : " << minVal << "     " << maxVal << std::endl;

    cv::normalize(image_dsm, image_dsm, 0, 1, cv::NORM_MINMAX, -1, image_dsm!=-1);
    cv::minMaxLoc(image_dsm, &minVal, &maxVal);
    std::cout << "Normalized image min-max values : " << minVal << "     " << maxVal << std::endl;

    cv::Mat image_dsm_resize;
    cv::resize(image_dsm, image_dsm_resize, cv::Size(40000, 39999), 0, 0, cv::INTER_CUBIC);
    cout << "Resized image size : " << image_dsm_resize.rows << " " << image_dsm_resize.cols << endl;
    cv::minMaxLoc(image_dsm_resize, &minVal, &maxVal);
    std::cout << "Resized image min-max values : " << minVal << "     " << maxVal << std::endl;

    std::string basepath = path.substr(0, path.find("."));
    for (int y=0; ((y+1)*1024+720)<=39999; y++) {
        for (int x=0; ((x+1)*1024)<=40000;x++) {
            cv::Range rows(y*1024+720, (y+1)*1024+720);
            cv::Range cols(x*1024, (x+1)*1024);
            cv::Mat image_dsm_cropped = image_dsm_resize(rows, cols);
            cv::imwrite(basepath+"_"+"CUB_"+to_string(y*1024)+"_"+to_string(x*1024)+".tif", image_dsm_cropped);
        }
    }

}
