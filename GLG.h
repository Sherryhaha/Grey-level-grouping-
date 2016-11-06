//
// Created by sunguoyan on 2016/11/4.
//

#ifndef GREYLG_GLG_H
#define GREYLG_GLG_H
#include <cv.h>
#include <iostream>
#include <opencv2/highgui.hpp>
#include <vector>
#include <opencv2/imgproc.hpp>
using namespace cv;
using namespace std;

typedef uchar tc;

class GLG{
public:
    int X_image,Y_image;
    int M=256;         //grayscale
    int n=0;         //nonzero histogram components

    void getHistogram(Mat&I,vector<int> &H);
    void getH(Mat&I,vector<int> &H);
    void reGroup(vector<vector<float>>&G,vector<vector<float>>&L,vector<vector<float>>&R,int gn);
    void transformation(int gn,vector<int>&T,vector<float>&L,vector<float>&R);
    void doMap(vector<int>&T,Mat&I,Mat&D);
    void GGLg(Mat&src,Mat&final);
    void fGLg(Mat&src,Mat&final);
    double distance(Mat&D);
    void show_Hist(Mat&src);
    double TEN(Mat&I,int T);
    double Tenegrad(Mat& src);
};


#endif //GREYLG_GLG_H
