//
// Created by sunguoyan on 2016/11/4.
//

#include "GLG.h"

//计算传递函数
void GLG::transformation(int cycle, vector<int> &T, vector<float> &L, vector<float> &R) {
    float N, alpha;
    if (L[0] == R[0])
        alpha = 0.8;
    else
        alpha = 0;
//    N = (M - 1) / (cycle - alpha);
    N = (M - 1) / (n - cycle - alpha);
    for (int k = 0; k < 256; k++) {
        if (k <= L[0]) {
            T[k] = 0;
//        } else if (k >= R[cycle-1]) {
        } else if (k >= R[n - cycle - 1]) {
            T[k] = 255;
        } else {
            for (int i = 0; i < n - cycle; i++) {
//                for (int i = 0; i < cycle; i++) {
                if ((k < R[i]) && (k >= L[i])) {
                    if (L[i] == R[i]) {
                        T[k] = (i + 1 - alpha) * N;
                    } else {
                        T[k] = (i + 1 - alpha - (R[i] - k) / (R[i] - L[i])) * N + 1;
                    }
                }
            }
            for (int i = 0; i < n - cycle + 1; i++) {
                if ((k < L[i + 1]) && (k >= R[i])) {
                    T[k] = (i + 1 - alpha) * N;
                }
            }
        }
    }
}


void GLG::GGLg(Mat &src, Mat &final) {
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

    //把G，L,R初始化
    int i = 0;
    for (int k = 0; k < 256; k++) {
        if (H[k] != 0) {
            G[0][i] = H[k];
            L[0][i] = k;
            R[0][i] = k;
            i++;
        }
    }


    int cycle = 20;
    vector<vector<int>> T(M);
    for (int i = 0; i < M; i++) {
        T[i].resize(M);
    }
    vector<vector<int>> Hn(M);
    for (int i = 0; i < M; i++) {
        Hn[i].resize(M);
    }
    for (int i = 0; i < M; i++) {
        Hn[0][i] = H[i];
    }

    for (int ii = 0; ii < n - cycle; ii++) {
        reGroup(G, L, R, ii);
    }
    double dis = 0, MD = 0;
    int c;
    //就算传递函数，并计算出最大距离
    for (int ii = 1; ii < n - cycle + 1; ii++) {
        transformation(ii, T[ii], L[ii], R[ii]);
        doMap(T[ii], src, final);
        dis = distance(final);
        if (dis > MD) {
            MD = dis;
            c = ii;
        }
    }
    cout << "nonzero histogram components: " << n << endl;
    cout << "diatance of the max processed image: " << MD << endl;
    double SD = distance(src);
    cout << "distance of the original image: " << SD << endl;
    doMap(T[c], src, final);
}
int main() {
    string filename = "/Users/sunguoyan/Downloads/picture/otherPicture/bag.png";
    Mat src, final, Hs, e;
    src = imread(filename, CV_LOAD_IMAGE_GRAYSCALE);
    Hs = src.clone();

    GLG t;
    t.X_image = src.cols;
    t.Y_image = src.rows;
    final.create(t.Y_image, t.X_image, CV_8UC1);

    t.GGLg(src, final);

    //直方图均衡化
    double D;
    equalizeHist(Hs, e);
    D = t.distance(e);
    cout << "distance of hisrogram: " << D << endl;
    t.show_Hist(e);
    double ten = t.Tenegrad(final);
    double sten = t.Tenegrad(src);
    double hten = t.Tenegrad(e);
    cout<<"ten of original image:"<<sten<<endl;

    cout<<"ten of HE processed image:"<<hten<<endl;
    cout<<"ten of gGLG processed image:"<<ten<<endl;

    namedWindow("original");
    namedWindow("afterProcess");
    imshow("original", src);
    imshow("afterProcess", e);
    waitKey(0);
    return 0;
}

