//
// Created by sunguoyan on 2016/11/4.
//

#include "GLG.h"

//得到图像的直方图，并算出非0值
void GLG::getHistogram(Mat &I, vector<int> &H) {
    for (int h = 0; h < Y_image; h++) {
        for (int w = 0; w < X_image; w++) {
            if (I.at<tc>(h, w) > M || I.at<tc>(h, w) < 0) {
                cout << "wrong" << endl;
            }
            H[I.at<tc>(h, w)]++;
        }
    }

    for (int i = 0; i < M; i++) {
        if (H[i] != 0) {
            n++;
        }
    }
}


//将图像的直方图进行重组,计算重组后的G,L,R
void GLG::reGroup(vector<vector<float>> &G, vector<vector<float>> &L, vector<vector<float>> &R, int r) {
    int a = G[r][0];
    int ia = 0;
    for (int i = 0; i < n - r; i++) {
        if (G[r][i] < a) {
            a = G[r][i];
            ia = i;
        }
    }
    int b;
    int j;

    if (ia == n - r - 1) {
        b = G[r][n - r - 2];
        j = n - r - 1;
    } else if (G[r][ia - 1] <= G[r][ia + 1]) {
        b = G[r][ia - 1];
        j = ia;
    } else {
        b = G[r][ia + 1];
        j = ia + 1;
    }

    //G计算
    for (int i = 0; i < j - 1; i++) {
        G[r + 1][i] = G[r][i];
    }
    G[r + 1][j - 1] = a + b;
    for (int i = j; i < n - r - 1; i++) {
        G[r + 1][i] = G[r][i + 1];
    }

    //L的计算
    for (int i = 0; i < j; i++) {
        L[r + 1][i] = L[r][i];
    }
    for (int i = j; i < n - r - 1; i++) {
        L[r + 1][i] = L[r][i + 1];
    }
    //R的计算
    for (int i = 0; i < j - 1; i++) {
        R[r + 1][i] = R[r][i];
    }
    for (int i = j - 1; i < n - r - 1; i++) {
        R[r + 1][i] = R[r][i + 1];
    }
}

////计算传递函数
//void GLG::transformation(int cycle, vector<int> &T, vector<float> &L, vector<float> &R) {
//    float N, alpha;
//    if (L[0] == R[0])
//        alpha = 0.8;
//    else
//        alpha = 0;
//    N = (M - 1) / (cycle - alpha);
//    for (int k = 0; k < 256; k++) {
//        if (k <= L[0]) {
//            T[k] = 0;
//        } else if (k >= R[cycle - 1]) {
//            T[k] = 255;
//        } else {
//            //在一个bins中
//            for (int i = 0; i < cycle; i++) {
//                if ((k <= R[i]) && (k >= L[i])) {
//                    //该bin只有一个灰度级
//                    if (L[i] == R[i]) {
//                        T[k] = (i + 1 - alpha) * N;
//                    }
//                        //该bin包含多个灰度级
//                    else {
//                        T[k] = (i + 1 - alpha - (R[i] - k) / (R[i] - L[i])) * N + 1;
//                    }
//                }
//            }
//            //在两个bins之间
//            for (int i = 0; i < cycle - 1; i++) {
//                if ((k < L[i + 1]) && (k > R[i])) {
//                    T[k] = (i + 1 - alpha) * N;
//                }
//            }
//        }
//    }
//}

//将图像按照传递函数进行映射
void GLG::doMap(vector<int> &T, Mat &I, Mat &D) {
    for (int h = 0; h < Y_image; h++) {
        for (int w = 0; w < X_image; w++) {
            D.at<tc>(h, w) = (uchar) (T[I.at<tc>(h, w)] + 0.5);
        }
    }
}

double GLG::distance(Mat &D) {
    vector<int> H(M, 0);
    const int channels[1] = {0};
    const int histSize[1] = {256};
    float hrange[2] = {0, 255};
    const float *range[1] = {hrange};

    MatND hist;
    cv::calcHist(&D, 1, channels, Mat(), hist, 1, histSize, range);
    long Dis = 0;
    for (int i = 0; i < M - 1; i++) {
        for (int j = i + 1; j < M; j++) {
            Dis += hist.at<float>(i) * hist.at<float>(j) * (j - i);
        }
    }
//    cout << Dis << endl;
    double Npix = X_image * Y_image;
    double N = (Npix) * ((Npix - 1));
//    cout << N << endl;
    double d;
    d = Dis / N;
    return d;
}

//绘制直方图
void GLG::show_Hist(Mat&src){
    int bins = 256;
    int hist_size[] = {bins};
    float range[] = { 0, 256 };
    const float* ranges[] = { range};
    MatND hist;
    int channels[] = {0};

    calcHist( &src, 1, channels, Mat(), // do not use mask
              hist, 1, hist_size, ranges,
              true, // the histogram is uniform
              false );

    double max_val;
    minMaxLoc(hist, 0, &max_val, 0, 0);
    int scale = 2;
    int hist_height=256;
    Mat hist_img = Mat::zeros(hist_height,bins*scale, CV_8UC3);
    for(int i=0;i<bins;i++)
    {
        float bin_val = hist.at<float>(i);
        int intensity = cvRound(bin_val*hist_height/max_val);  //要绘制的高度
        rectangle(hist_img,Point(i*scale,hist_height-1),
                  Point((i+1)*scale - 1, hist_height - intensity),
                  CV_RGB(255,255,255));
    }
    imshow( "Histogram proccessed", hist_img );
}

double GLG::TEN(Mat&I,int T){
    Mat image=I.clone();

    double tmpx,tmpy;
    for(int i = 1;i<Y_image-1;i++){
        for(int j = 1;j<X_image-1;j++){
            tmpx=abs(I.at<tc>(i-1,j+1)+I.at<tc>(i,j+1)*2+I.at<tc>(i+1,j+1)-I.at<tc>(i-1,j-1)-I.at<tc>(i,j-1)*2-I.at<tc>(i+1,j-1));
            tmpy = abs(I.at<tc>(i-1,j-1)+I.at<tc>(i-1,j)*2+I.at<tc>(i-1,j+1)-I.at<tc>(i+1,j-1)-I.at<tc>(i+1,j)*2-I.at<tc>(i+1,j+1));
//            image.at<tc>(i,j) = sqrt(pow(tmpx*image.at<tc>(i,j),2)+pow(tmpy*image.at<tc>(i,j),2));
            image.at<tc>(i,j) = sqrt(pow(tmpx,2)+pow(tmpy,2));

        }
    }
    double tmp=0;
    for(int i = 1;i<Y_image-1;i++){
        for(int j = 1;j<X_image-1;j++){
            if(image.at<tc>(i,j)>T){
                tmp+=image.at<tc>(i,j);
            }
        }
    }
    return tmp;
}


double GLG::Tenegrad(Mat& src)
{
    int widthstep=src.step;
    uchar *data=src.data;
    double S=0;
    for(int x = 1;x<Y_image-1;x++)
    {
        uchar *pre_row=data +(x-1)*widthstep;
        uchar *cur_row=data +x*widthstep;
        uchar *nex_row=data +(x+1)*widthstep;
        int Sx,Sy;
        for(int y = 1;y<X_image-1;y++)
        {
            Sx=(uchar)pre_row[y+1]+2*(uchar)cur_row[y+1]+(uchar)nex_row[y+1]
               -(uchar)pre_row[y-1]-2*(uchar)cur_row[y-1]-(uchar)nex_row[y-1];
            Sy=(uchar)nex_row[y-1]+2*(uchar)nex_row[y]+(uchar)nex_row[y+1]
               -(uchar)pre_row[y-1]-2*(uchar)pre_row[y]-(uchar)pre_row[y+1];
            S+=Sx*Sx+Sy*Sy;
        }
    }
    return S;
}