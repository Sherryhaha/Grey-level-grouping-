
//
// Created by sunguoyan on 2016/11/3.
//

#include "GLG.h"


//计算传递函数
void GLG::transformation(int cycle, vector<int> &T, vector<float> &L, vector<float> &R) {
    float N, alpha;
    if (L[0] == R[0])
        alpha = 0.8;
    else
        alpha = 0;
    N = (M - 1) / (cycle - alpha);
    for (int k = 0; k < 256; k++) {
        if (k <= L[0]) {
            T[k] = 0;
        } else if (k >= R[cycle - 1]) {
            T[k] = 255;
        } else {
            //在一个bins中
            for (int i = 0; i < cycle; i++) {
                if ((k <= R[i]) && (k >= L[i])) {
                    //该bin只有一个灰度级
                    if (L[i] == R[i]) {
                        T[k] = (i + 1 - alpha) * N;
                    }
                        //该bin包含多个灰度级
                    else {
                        T[k] = (i + 1 - alpha - (R[i] - k) / (R[i] - L[i])) * N + 1;
                    }
                }
            }
            //在两个bins之间
            for (int i = 0; i < cycle - 1; i++) {
                if ((k < L[i + 1]) && (k > R[i])) {
                    T[k] = (i + 1 - alpha) * N;
                }
            }
        }
    }
}



void GLG::fGLg(Mat &src, Mat &final) {
    vector<int> H(M);
    getHistogram(src, H);

    vector<vector<float>> G(M);
    vector<vector<float>> L(M);
    vector<vector<float>> R(M);
    for (int i = 0; i < M; i++) {
        G[i].resize(M);
        L[i].resize(M);
        R[i].resize(M);
    }
    int i = 0;
    //把G，L,R初始化
    for (int k = 0; k < 256; k++) {
        if (H[k] != 0) {
            G[0][i] = H[k];
            L[0][i] = k;
            R[0][i] = k;
            i++;
        }
    }

    int cycle = 20;

    vector<int> TT(M, 0);

    for (int ii = 0; ii < n - cycle; ii++) {
        reGroup(G, L, R, ii);
    }

    transformation(cycle, TT, L[n - cycle], R[n - cycle]);

    doMap(TT, src, final);
    double d, sd;
    sd = distance(src);
    d = distance(final);
    cout << "distance of the original image: " << sd << endl;
    cout << "distance of the processed image: " << d << endl;
}



int main() {
    string filename = "/Users/sunguoyan/Downloads/picture/otherPicture/bag.png";
    Mat src, final;
    src = imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
    GLG t;
    t.X_image = src.cols;
    t.Y_image = src.rows;
    final.create(t.Y_image, t.X_image, CV_8UC1);


    t.fGLg(src, final);


    namedWindow("original");
    namedWindow("afterProcess");
    imshow("original", src);
    imshow("afterProcess", final);
    double ften = t.Tenegrad(final);
    cout<<"ten of fast GLG image:"<<ften<<endl;

    MatND hist2;
//    myCal_Hist(final, hist2);
    t.show_Hist(final);

    waitKey(0);

    return 0;
}