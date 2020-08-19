#include "slic6d.h"

Slic6D::Slic6D():num_superpixels(0), step(0)
{

}

Slic6D::~Slic6D()
{

}

void Slic6D::setDsm(cv::Mat *value)
{
    dsm = value;
}

cv::Mat *Slic6D::getConn_clusters() const
{
    return conn_clusters;
}

void Slic6D::init()
{
    //std::cout << "Start : Init method." << std::endl;

    clusters = new cv::Mat(image->rows, image->cols, CV_16SC1);
    clusters->setTo(-1);
    //std::cout << "Cluster matrix size : " << clusters->rows << " " << clusters->cols << std::endl;
    //std::cout << "Cluster matrix init value : " << clusters->at<short>(1000,1000) << std::endl;

    distances = new cv::Mat(image->rows, image->cols, CV_32FC1);
    distances->setTo(DBL_MAX);
    //std::cout << "Distance matrix size : " << distances->rows << " " << distances->cols << std::endl;
    //std::cout << "Distance matrix init value : " << distances->at<float>(67,67) << std::endl;

    cv::Mat *gray_image = new cv::Mat;
    cv::cvtColor(*image, *gray_image, cv::COLOR_BGR2GRAY);

    step = sqrt ( (image->rows * image->cols) / num_superpixels );
    for (int y = step/2; y < image->rows - step/2; y += step) {
        for (int x = step/2; x < image->cols - step/2; x += step) {
            std::vector<float> center;

            /* Find the local minimum (gradient-wise). */
            cv::Point2i min_center = find_local_minimum(gray_image, cv::Point2i(x,y));
            float min_center_dsm = dsm->at<float>(min_center.y, min_center.x);
            cv::Vec3b min_center_colour = image->at<cv::Vec3b>(min_center.y, min_center.x);

            /* Generate the center vector. */
            center.push_back((float)min_center_colour.val[0]);
            center.push_back((float)min_center_colour.val[1]);
            center.push_back((float)min_center_colour.val[2]);
            center.push_back((float)min_center.x);
            center.push_back((float)min_center.y);
            center.push_back(min_center_dsm);
            center.push_back(0.0);

            /* Append to vector of centers. */
            centers.push_back(center);
            //center_counts.push_back(0);
        }
    }
    //std::cout << "End : Init method" << std::endl;
}

void Slic6D::generate_superpixels(cv::Mat *input_image, cv::Mat *dsm, int num_superpixels)
{
    setImage(input_image);
    setDsm(dsm);
    setNum_superpixels(num_superpixels);
    //std::cout << "Image size : " << input_image->rows << " " << input_image->cols << std::endl;

    init();

    cv::Vec3b color;
    cv::Point3f coordinate;
    float distance;
    for (int iter=0; iter<25; iter++) {

        distances->setTo(DBL_MAX);

        for (short c_index=0; c_index<centers.size(); c_index++) {
            for (int y=centers[c_index][4]-step; y<=centers[c_index][4]+step; y++) {
                for (int x=centers[c_index][3]-step; x<=centers[c_index][3]+step; x++) {
                    if(y>=0 && y<image->rows && x>=0 && x<image->cols){
                        color = image->at<cv::Vec3b>(y,x);

                        coordinate.x = x;
                        coordinate.y = y;
                        coordinate.z = dsm->at<float>(y,x);

                        distance = compute_dist(c_index, color, coordinate);

                        if(distance < distances->at<float>(y,x)){
                            distances->at<float>(y,x) = distance;
                            clusters->at<short>(y,x) = c_index;
                        }
                    }
                }
            }
        }

        /* Clear the center values. */
        for (int i = 0; i < (int) centers.size(); i++) {
            centers[i][0] = centers[i][1] = centers[i][2] = centers[i][3] = centers[i][4] = centers[i][5] = centers[i][6] = 0;
        }

        /* Compute the new cluster centers. */
        for (int y=0; y<image->rows; y++) {
            for (int x=0; x<image->cols; x++) {
                short c_id = clusters->at<short>(y,x);

                if (c_id != -1) {
                    cv::Vec3b colour = image->at<cv::Vec3b>(y,x);

                    centers[c_id][0] += (float)colour.val[0];
                    centers[c_id][1] += (float)colour.val[1];
                    centers[c_id][2] += (float)colour.val[2];
                    centers[c_id][3] += x;
                    centers[c_id][4] += y;
                    centers[c_id][5] += this->dsm->at<float>(y,x);
                    centers[c_id][6] += 1;

                }
            }
        }

        /* Normalize the clusters. */
        for (int i = 0; i < centers.size(); i++) {
            int center_size = centers[i].size();
            for (int j = 0; j < center_size-1; j++){
                centers[i][j] /= centers[i][center_size-1];
            }
        }

    }

//    cv::Mat clustercon;
//    clusters->convertTo(clustercon, CV_8U);
//    clusters->setTo(0, *clusters!=15);
//    clusters->setTo(50, *clusters==98);
//    clusters->setTo(50, *clusters==97);
//    clusters->setTo(50, *clusters==96);
//    double minVal, maxVal;
//    cv::minMaxLoc(*clusters, &minVal, &maxVal);
//    std::cout << "Cluster min-max values : " << minVal << "     " << maxVal << std::endl;
 //   cv::normalize(*clusters, clustercon, 0, 32767, cv::NORM_MINMAX);
//    cv::imshow("Clusters", clustercon);
//    cv::waitKey();
 //   display_clusters("Initial cluster", false);
    //std::cout << "Number of clusters : " << centers.size() << std::endl;
}

void Slic6D::create_connectivity()
{
    int min_cluster_size = (image->rows*image->cols)/(this->num_superpixels*10);
    conn_clusters = new cv::Mat(image->rows, image->cols, CV_16SC1);
    conn_clusters->setTo(-1);

    int dy4[4] = { 0,  1,  0, -1};
    int dx4[4] = { 1,  0, -1,  0};

    short label=0, adj_label=0;

    for (int y=0; y<image->rows; y++) {
        for (int x=0; x<image->cols; x++) {
            if(conn_clusters->at<short>(y,x)==-1){
                std::vector<cv::Point2i> cluster_points;
                cluster_points.push_back(cv::Point(x,y));

                /* Find an adjacent label, for possible use later. */
                for (int d=0; d<4; d++) {
                    int dx = x + dx4[d], dy = y + dy4[d];
                    if (dx >= 0 && dx < image->cols && dy >= 0 && dy < image->rows) {
                        if (conn_clusters->at<short>(dy,dx) >= 0) {
                            adj_label = conn_clusters->at<short>(dy,dx);
                        }
                    }
                }

                for (int cp_index=0; cp_index<cluster_points.size();cp_index++) {
                    for (int d=0; d<4; d++) {
                        int dy = cluster_points[cp_index].y+dy4[d], dx = cluster_points[cp_index].x+dx4[d];
                        if(dy>=0 && dy<image->rows && dx>=0 && dx<image->cols){
                            if(conn_clusters->at<short>(dy,dx)==-1 &&
                                    (clusters->at<short>(dy,dx)==clusters->at<short>(y,x)) ){
                                conn_clusters->at<short>(dy,dx) = label;
                                cluster_points.push_back(cv::Point2i(dx,dy));
                            }
                        }
                    }
                }
                label += 1;

                if(cluster_points.size()<min_cluster_size){
                    for (int cp_index=0;cp_index<cluster_points.size();cp_index++) {
                        conn_clusters->at<short>(cluster_points[cp_index].y,cluster_points[cp_index].x) = adj_label;
                    }
                    label -=1;
                }
            }
        }
    }
    //std::cout << "Number of connected clusters : " << label << std::endl;
}


cv::Mat Slic6D::display_clusters(std::string win_name, bool connected)
{
    cv::Mat *cluster;
    if(connected==true){
        cluster = conn_clusters;
    }
    else {
        cluster = clusters;
    }

    cv::Mat disp_cluster = cv::Mat::zeros(cv::Size(image->cols,image->rows), CV_8UC3);
    //short r=0,g=0,b=0;
    for (int y=0; y<this->image->rows; y++) {
        for (int x=0; x<this->image->cols; x++) {
            short val = cluster->at<short>(y,x);
            if(val%3 == 0){
                //r = (r+30)%256;
                disp_cluster.at<cv::Vec3b>(y,x) = cv::Vec3b(255,0,0);
            }
            else if (val%3 == 1) {
                //g = (g+30)%256;
                disp_cluster.at<cv::Vec3b>(y,x) = cv::Vec3b(0,255,0);
            }
            else if (val%3 == 2) {
                //b = (b+30)%256;
                disp_cluster.at<cv::Vec3b>(y,x) = cv::Vec3b(0,0,255);
            }
        }
    }
    cv::Mat disp_image ;//= cv::Mat::zeros(cv::Size(image->cols,image->rows), CV_8SC3);
    //std::cout << image->type() << " " << disp_cluster.type() << std::endl;


    cv::Mat edge;
    cluster->convertTo(edge,CV_8U);
    cv::Canny(edge, edge, 0, 0*3, 3);
    cv::dilate(edge, edge, cv::getStructuringElement(cv::MORPH_RECT,cv::Size(3,3)));
//    cv::imshow("Edge", edge);
//    cv::waitKey();

    disp_image = image->clone();
    //cv::addWeighted(*image, 0.90, disp_cluster, 0.10, 0, disp_image, CV_8UC3);
    disp_image.setTo(cv::Vec3b(0,0,255), edge==255);
    //cv::imshow(win_name, disp_image);
    //cv::waitKey();
    return disp_image;
}

float Slic6D::compute_dist(int center_index, cv::Vec3b color, cv::Point3f coordinate)
{
    float alpha_color = 3.0/255.0;
    float alpha_coord = 1.0/step;
    float alpha_dsm = 1.0;
    float dist_color = sqrt(pow(centers[center_index][0] - color.val[0], 2)
                           + pow(centers[center_index][1] - color.val[1], 2)
                           + pow(centers[center_index][2] - color.val[2], 2));
    float dist_coord = sqrt(pow(centers[center_index][3] - coordinate.x, 2)
                           + pow(centers[center_index][4] - coordinate.y, 2));
                           //+ pow(centers[center_index][5] - coordinate.z, 2));
    float dist_dsm = sqrt(pow(centers[center_index][5] - coordinate.z, 2));

    //float distance = dist_color + dist_coord;
    float distance = alpha_color*dist_color +
                     alpha_coord*dist_coord +
                     alpha_dsm*dist_dsm;
    return distance;
    //return sqrt(pow(dist_color / 1.0, 2) + pow(dist_coord / 1.0, 2));
    //return 0.0;
}

cv::Point2i Slic6D::find_local_minimum(cv::Mat *image, cv::Point2i center)
{
    float min_grad = FLT_MAX;
    cv::Point2i loc_min = center;

    for (int y = center.y-1; y < center.y+2; y++) {
        for (int x = center.x-1; x < center.x+2; x++) {
            float c_pixel   = (float)image->at<uchar>(y,x);
            float c_pixel_y = (float)image->at<uchar>(y+1,x);
            float c_pixel_x = (float)image->at<uchar>(y,x+1);

            /* Compute horizontal and vertical gradients and keep track of the
               minimum. */
            if((abs(c_pixel_y - c_pixel) + abs(c_pixel_x - c_pixel)) < min_grad){
                min_grad = abs(c_pixel_y - c_pixel) + abs(c_pixel_x - c_pixel);
                loc_min.x = x;
                loc_min.y = y;
            }
        }
    }

    return loc_min;
}

void Slic6D::setNum_superpixels(int value)
{
    num_superpixels = value;
}

void Slic6D::setImage(cv::Mat *value)
{
    image = value;
}
