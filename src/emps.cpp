#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define IN_FILE "../mk_particle/dambreak.prof"
#define PCL_DST 0.02                    //平均粒子間距離
#define MIN_X  (0.0 - PCL_DST*3)    //解析領域のx方向の最小値
#define MIN_Y  (0.0 - PCL_DST*3)    //解析領域のy方向の最小値
#define MIN_Z  (0.0 - PCL_DST*3)    //解析領域のz方向の最小値
#define MAX_X  (1.0 + PCL_DST*3)    //解析領域のx方向の最大値
#define MAX_Y  (0.2 + PCL_DST*3)    //解析領域のy方向の最大値
#define MAX_Z  (0.6 + PCL_DST*30)    //解析領域のz方向の最大値

#define GST -1            //計算対象外粒子の種類番号
#define FLD 0                //流体粒子の種類番号
#define WLL  1            //壁粒子の種類番号
#define NUM_TYP  2        //粒子の種類数
#define DNS_FLD 1000        //流体粒子の密度
#define DNS_WLL 1000        //壁粒子の密度
#define DT 0.0005            //時間刻み幅
#define FIN_TIM 1.0        //時間の上限
#define SND 22.0            //音速
#define OPT_FQC 100        //出力間隔を決める反復数
#define KNM_VSC 0.000001    //動粘性係数
#define DIM 3                //次元数
#define CRT_NUM 0.1        //クーラン条件数
#define COL_RAT 0.2        //接近した粒子の反発率
#define DST_LMT_RAT 0.9    //これ以上の粒子間の接近を許さない距離の係数
#define G_X 0.0            //重力加速度のx成分
#define G_Y 0.0            //重力加速度のy成分
#define G_Z -9.8            //重力加速度のz成分
#define WEI(dist, re) ((re/dist) - 1.0)    //重み関数

FILE *fp;
char filename[256];
int iLP, iF;
double TIM;
int nP;
double *Acc, *Pos, *Vel, *Prs, *pav;
int *Typ;
double r, r2;
double DB, DB2, DBinv;
int nBx, nBy, nBz, nBxy, nBxyz;
int *bfst, *blst, *nxt;
double n0, lmd, A1, A2, A3, rlim, rlim2, COL;
double Dns[NUM_TYP], invDns[NUM_TYP];

//関数01 計算領域を飛び出たら幽霊粒子に設定し計算から排除
void ChkPcl(int i) {
    if (Pos[i * 3] > MAX_X || Pos[i * 3] < MIN_X ||
        Pos[i * 3 + 1] > MAX_Y || Pos[i * 3 + 1] < MIN_Y ||
        Pos[i * 3 + 2] > MAX_Z || Pos[i * 3 + 2] < MIN_Z) {
        Typ[i] = GST;
        Prs[i] = Vel[i * 3] = Vel[i * 3 + 1] = Vel[i * 3 + 2] = 0.0;
    }
}

//関数02 初期状態読み込み
void RdDat(void) {
    fp = fopen(IN_FILE, "r");
    fscanf(fp, "%d", &nP);
    printf("nP: %d\n", nP);
    Acc = (double *) malloc(sizeof(double) * nP * 3);    //粒子の加速度
    Pos = (double *) malloc(sizeof(double) * nP * 3);    //粒子の座標
    Vel = (double *) malloc(sizeof(double) * nP * 3);    //粒子の速度
    Prs = (double *) malloc(sizeof(double) * nP);        //粒子の圧力
    pav = (double *) malloc(sizeof(double) * nP);        //時間平均された粒子の圧力
    Typ = (int *) malloc(sizeof(int) * nP);            //粒子の種類
    for (int i = 0; i < nP; i++) {
        int a[2];
        double b[8];
        //粒子番号　粒子種別(GST -1 FLD 0 WLL) 座標三次元　速度三次元　（加速度は0）
        fscanf(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf", &a[0], &a[1], &b[0], &b[1], &b[2], &b[3], &b[4], &b[5],
               &b[6], &b[7]);
        Typ[i] = a[1];
        Pos[i * 3] = b[0];
        Pos[i * 3 + 1] = b[1];
        Pos[i * 3 + 2] = b[2];
        Vel[i * 3] = b[3];
        Vel[i * 3 + 1] = b[4];
        Vel[i * 3 + 2] = b[5];
        Prs[i] = b[6];
        pav[i] = b[7];
    }
    fclose(fp);
    for (int i = 0; i < nP; i++) { ChkPcl(i); }
    for (int i = 0; i < nP * 3; i++) { Acc[i] = 0.0; }
}

//関数03 状態書き込み
void WrtDat(void) {
    char outout_filename[256];
    sprintf(outout_filename, "output%05d.prof", iF);
    fp = fopen(outout_filename, "w");
    fprintf(fp, "%d\n", nP);
    for (int i = 0; i < nP; i++) {
        int a[2];
        double b[8];
        a[0] = i;
        a[1] = Typ[i];
        b[0] = Pos[i * 3];
        b[1] = Pos[i * 3 + 1];
        b[2] = Pos[i * 3 + 2];
        b[3] = Vel[i * 3];
        b[4] = Vel[i * 3 + 1];
        b[5] = Vel[i * 3 + 2];
        b[6] = Prs[i];
        b[7] = pav[i] / OPT_FQC;
        fprintf(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", a[0], a[1], b[0], b[1], b[2], b[3], b[4], b[5], b[6],
                b[7]);
        //書き込みのタイミングで時間平均圧力をリセットしている
        pav[i] = 0.0;
    }
    fclose(fp);
    iF++;
}

//関数04 バケット構造のメモリを確保する
void AlcBkt(void) {
    //なんで2.1なの？
    r = PCL_DST * 2.1;        //影響半径
    r2 = r * r;
    DB = r * (1.0 + CRT_NUM);    //バケット1辺の長さ
    DB2 = DB * DB;
    DBinv = 1.0 / DB;
    nBx = (int) ((MAX_X - MIN_X) * DBinv) + 3;//解析領域内のx方向のバケット数
    nBy = (int) ((MAX_Y - MIN_Y) * DBinv) + 3;//解析領域内のy方向のバケット数
    nBz = (int) ((MAX_Z - MIN_Z) * DBinv) + 3;//解析領域内のz方向のバケット数
    nBxy = nBx * nBy;
    nBxyz = nBx * nBy * nBz;        //解析領域内のバケット数
    printf("nBx:%d  nBy:%d  nBz:%d  nBxy:%d  nBxyz:%d\n", nBx, nBy, nBz, nBxy, nBxyz);
    bfst = (int *) malloc(sizeof(int) * nBxyz);    //バケットに格納された先頭の粒子番号
    blst = (int *) malloc(sizeof(int) * nBxyz);    //バケットに格納された最後尾の粒子番号
    nxt = (int *) malloc(sizeof(int) * nP);        //同じバケット内の次の粒子番号
}

//関数05 粒子法にかかわるパラメータの設定
void SetPara(void) {
    //理想的な格子の状態で粒子数密度とラプラシアンモデルの係数を求めている
    double tn0 = 0.0;
    double tlmd = 0.0;
    for (int ix = -4; ix < 5; ix++) {
        for (int iy = -4; iy < 5; iy++) {
            for (int iz = -4; iz < 5; iz++) {
                double x = PCL_DST * (double) ix;
                double y = PCL_DST * (double) iy;
                double z = PCL_DST * (double) iz;
                double dist2 = x * x + y * y + z * z;
                if (dist2 <= r2) {
                    if (dist2 == 0.0)continue;
                    double dist = sqrt(dist2);
                    tn0 += WEI(dist, r);
                    tlmd += dist2 * WEI(dist, r);
                }
            }
        }
    }
    n0 = tn0;            //初期粒子数密度
    lmd = tlmd / tn0;    //ラプラシアンモデルの係数λ

    A1 = 2.0 * KNM_VSC * DIM / n0 / lmd;//粘性項の計算に用いる係数
    A2 = SND * SND / n0;                //圧力の計算に用いる係数
    A3 = -DIM / n0;                    //圧力勾配項の計算に用いる係数
    Dns[FLD] = DNS_FLD;
    Dns[WLL] = DNS_WLL;
    invDns[FLD] = 1.0 / DNS_FLD;
    invDns[WLL] = 1.0 / DNS_WLL;
    rlim = PCL_DST * DST_LMT_RAT;//これ以上の粒子間の接近を許さない距離
    rlim2 = rlim * rlim;
    COL = 1.0 + COL_RAT;
    iLP = 0;            //反復数
    iF = 0;            //ファイル番号
    TIM = 0.0;        //時刻
}

//関数06 バケットのxyzを設定する
void MkBkt(void) {
    for (int i = 0; i < nBxyz; i++) { bfst[i] = -1; }
    for (int i = 0; i < nBxyz; i++) { blst[i] = -1; }
    for (int i = 0; i < nP; i++) { nxt[i] = -1; }
    for (int i = 0; i < nP; i++) {
        //計算しない
        if (Typ[i] == GST)continue;
        //DBinv = 1.0/DB;
        int ix = (int) ((Pos[i * 3] - MIN_X) * DBinv) + 1;
        int iy = (int) ((Pos[i * 3 + 1] - MIN_Y) * DBinv) + 1;
        int iz = (int) ((Pos[i * 3 + 2] - MIN_Z) * DBinv) + 1;
        //三次元的番号ib
        int ib = iz * nBxy + iy * nBx + ix;
        int j = blst[ib];
        //メモリ確保を最小限にするために特殊なリストを使っている、bfstは同バケット内の最初の粒子のリスト、blasは同バケット内の最後の粒子のリスト、nxtは同バケット内の次粒子のリスト
        blst[ib] = i;
        if (j == -1) { bfst[ib] = i; }
        else { nxt[j] = i; }
    }
}

//関数07 粘性項と重力を計算する関数
void VscTrm() {
    for (int i = 0; i < nP; i++) {
        if (Typ[i] == FLD) {
            double Acc_x = 0.0;
            double Acc_y = 0.0;
            double Acc_z = 0.0;
            double pos_ix = Pos[i * 3];
            double pos_iy = Pos[i * 3 + 1];
            double pos_iz = Pos[i * 3 + 2];
            double vec_ix = Vel[i * 3];
            double vec_iy = Vel[i * 3 + 1];
            double vec_iz = Vel[i * 3 + 2];
            int ix = (int) ((pos_ix - MIN_X) * DBinv) + 1;
            int iy = (int) ((pos_iy - MIN_Y) * DBinv) + 1;
            int iz = (int) ((pos_iz - MIN_Z) * DBinv) + 1;
            for (int jz = iz - 1; jz <= iz + 1; jz++) {
                for (int jy = iy - 1; jy <= iy + 1; jy++) {
                    for (int jx = ix - 1; jx <= ix + 1; jx++) {
                        int jb = jz * nBxy + jy * nBx + jx;
                        int j = bfst[jb];
                        if (j == -1) continue;
                        for (;;) {
                            //ここで周辺粒子が計算される
                            double v0 = Pos[j * 3] - pos_ix;
                            double v1 = Pos[j * 3 + 1] - pos_iy;
                            double v2 = Pos[j * 3 + 2] - pos_iz;
                            double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                            if (dist2 < r2) {
                                //近傍
                                if (j != i && Typ[j] != GST) {
                                    //自分以外の流体粒子
                                    //ウェイトと掛け算してAccに加算、Accの次元はまだ速度であるはず
                                    double dist = sqrt(dist2);
                                    double w = WEI(dist, r);
                                    Acc_x += (Vel[j * 3] - vec_ix) * w;
                                    Acc_y += (Vel[j * 3 + 1] - vec_iy) * w;
                                    Acc_z += (Vel[j * 3 + 2] - vec_iz) * w;
                                }
                            }
                            //イテレーター
                            j = nxt[j];
                            if (j == -1) break;
                        }
                    }
                }
            }
            //係数A1によってAccの次元が速度から加速度に変わる、また重力を加算
            Acc[i * 3] = Acc_x * A1 + G_X;
            Acc[i * 3 + 1] = Acc_y * A1 + G_Y;
            Acc[i * 3 + 2] = Acc_z * A1 + G_Z;
        }
    }
}

//関数08 加速度の修正項から位置と速度を求める関数
void UpPcl1() {
    for (int i = 0; i < nP; i++) {
        if (Typ[i] == FLD) {
            Vel[i * 3] += Acc[i * 3] * DT;
            Vel[i * 3 + 1] += Acc[i * 3 + 1] * DT;
            Vel[i * 3 + 2] += Acc[i * 3 + 2] * DT;
            Pos[i * 3] += Vel[i * 3] * DT;
            Pos[i * 3 + 1] += Vel[i * 3 + 1] * DT;
            Pos[i * 3 + 2] += Vel[i * 3 + 2] * DT;
            Acc[i * 3] = Acc[i * 3 + 1] = Acc[i * 3 + 2] = 0.0;
            //計算領域を飛び出たら幽霊粒子に設定し計算から排除
            ChkPcl(i);
        }
    }
}

//関数09 粒子の剛体衝突判定
void ChkCol() {
    for (int i = 0; i < nP; i++) {
        if (Typ[i] == FLD) {
            double mi = Dns[Typ[i]];
            double pos_ix = Pos[i * 3];
            double pos_iy = Pos[i * 3 + 1];
            double pos_iz = Pos[i * 3 + 2];
            double vec_ix = Vel[i * 3];
            double vec_iy = Vel[i * 3 + 1];
            double vec_iz = Vel[i * 3 + 2];
            double vec_ix2 = Vel[i * 3];
            double vec_iy2 = Vel[i * 3 + 1];
            double vec_iz2 = Vel[i * 3 + 2];
            int ix = (int) ((pos_ix - MIN_X) * DBinv) + 1;
            int iy = (int) ((pos_iy - MIN_Y) * DBinv) + 1;
            int iz = (int) ((pos_iz - MIN_Z) * DBinv) + 1;
            for (int jz = iz - 1; jz <= iz + 1; jz++) {
                for (int jy = iy - 1; jy <= iy + 1; jy++) {
                    for (int jx = ix - 1; jx <= ix + 1; jx++) {
                        int jb = jz * nBxy + jy * nBx + jx;
                        int j = bfst[jb];
                        if (j == -1) continue;
                        //周辺バケットの例のイテレータ
                        for (;;) {
                            //相対座標v=xj-xi
                            double v0 = Pos[j * 3] - pos_ix;
                            double v1 = Pos[j * 3 + 1] - pos_iy;
                            double v2 = Pos[j * 3 + 2] - pos_iz;
                            double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                            if (dist2 < rlim2) {
                                //衝突判定用二乗距離rlim2以内
                                if (j != i && Typ[j] != GST) {
                                    //かつ自分以外の計算可能粒子
                                    double fDT = (vec_ix - Vel[j * 3]) * v0 + (vec_iy - Vel[j * 3 + 1]) * v1 +
                                                 (vec_iz - Vel[j * 3 + 2]) * v2;
                                    if (fDT > 0.0) {
                                        //veli-velj * xj-xi >0 で内積から対向を判定
                                        //粒子種別に応じた密度
                                        double mj = Dns[Typ[j]];
                                        // TODO: この反発式は何？
                                        fDT *= COL * mj / (mi + mj) / dist2;
                                        vec_ix2 -= v0 * fDT;
                                        vec_iy2 -= v1 * fDT;
                                        vec_iz2 -= v2 * fDT;
                                    }
                                }
                            }
                            j = nxt[j];
                            if (j == -1) break;
                        }
                    }
                }
            }
            Acc[i * 3] = vec_ix2;
            Acc[i * 3 + 1] = vec_iy2;
            Acc[i * 3 + 2] = vec_iz2;
        }
    }
    for (int i = 0; i < nP; i++) {
        Vel[i * 3] = Acc[i * 3];
        Vel[i * 3 + 1] = Acc[i * 3 + 1];
        Vel[i * 3 + 2] = Acc[i * 3 + 2];
    }
}

//関数10 仮の圧力を求める
void MkPrs() {
    for (int i = 0; i < nP; i++) {
        if (Typ[i] != GST) {
            double pos_ix = Pos[i * 3];
            double pos_iy = Pos[i * 3 + 1];
            double pos_iz = Pos[i * 3 + 2];
            double ni = 0.0;
            int ix = (int) ((pos_ix - MIN_X) * DBinv) + 1;
            int iy = (int) ((pos_iy - MIN_Y) * DBinv) + 1;
            int iz = (int) ((pos_iz - MIN_Z) * DBinv) + 1;
            for (int jz = iz - 1; jz <= iz + 1; jz++) {
                for (int jy = iy - 1; jy <= iy + 1; jy++) {
                    for (int jx = ix - 1; jx <= ix + 1; jx++) {
                        int jb = jz * nBxy + jy * nBx + jx;
                        int j = bfst[jb];
                        if (j == -1) continue;
                        for (;;) {
                            //例のリスト構造、自分以外の計算可能な流体粒子
                            double v0 = Pos[j * 3] - pos_ix;
                            double v1 = Pos[j * 3 + 1] - pos_iy;
                            double v2 = Pos[j * 3 + 2] - pos_iz;
                            double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                            if (dist2 < r2) {
                                if (j != i && Typ[j] != GST) {
                                    double dist = sqrt(dist2);
                                    double w = WEI(dist, r);
                                    //単純に1*wを足している、自分を計算から省いているが圧力を求めるうえで適切なのか？
                                    ni += w;
                                }
                            }
                            j = nxt[j];
                            if (j == -1) break;
                        }
                    }
                }
            }
            double mi = Dns[Typ[i]];
            //技巧的な式、0か仮圧力が求まる、0は表面粒子、粒子毎にprs[i]に圧力を入れる
            double pressure = (ni > n0) * (ni - n0) * A2 * mi;
            Prs[i] = pressure;
        }
    }
}

//関数11 圧力勾配項を求める関数
void PrsGrdTrm() {
    for (int i = 0; i < nP; i++) {
        if (Typ[i] == FLD) {
            double Acc_x = 0.0;
            double Acc_y = 0.0;
            double Acc_z = 0.0;
            double pos_ix = Pos[i * 3];
            double pos_iy = Pos[i * 3 + 1];
            double pos_iz = Pos[i * 3 + 2];
            double pre_min = Prs[i];
            int ix = (int) ((pos_ix - MIN_X) * DBinv) + 1;
            int iy = (int) ((pos_iy - MIN_Y) * DBinv) + 1;
            int iz = (int) ((pos_iz - MIN_Z) * DBinv) + 1;
            for (int jz = iz - 1; jz <= iz + 1; jz++) {
                for (int jy = iy - 1; jy <= iy + 1; jy++) {
                    for (int jx = ix - 1; jx <= ix + 1; jx++) {
                        int jb = jz * nBxy + jy * nBx + jx;
                        int j = bfst[jb];
                        if (j == -1) continue;
                        for (;;) {
                            double v0 = Pos[j * 3] - pos_ix;
                            double v1 = Pos[j * 3 + 1] - pos_iy;
                            double v2 = Pos[j * 3 + 2] - pos_iz;
                            double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                            if (dist2 < r2) {
                                if (j != i && Typ[j] != GST) {
                                    //自分以外の周辺粒子の影響範囲計算可能粒子に最低圧力を保証している
                                    if (pre_min > Prs[j])pre_min = Prs[j];
                                }
                            }
                            j = nxt[j];
                            if (j == -1) break;
                        }
                    }
                }
            }
            for (int jz = iz - 1; jz <= iz + 1; jz++) {
                for (int jy = iy - 1; jy <= iy + 1; jy++) {
                    for (int jx = ix - 1; jx <= ix + 1; jx++) {
                        int jb = jz * nBxy + jy * nBx + jx;
                        int j = bfst[jb];
                        if (j == -1) continue;
                        for (;;) {
                            double v0 = Pos[j * 3] - pos_ix;
                            double v1 = Pos[j * 3 + 1] - pos_iy;
                            double v2 = Pos[j * 3 + 2] - pos_iz;
                            double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                            if (dist2 < r2) {
                                if (j != i && Typ[j] != GST) {
                                    //自分以外の周辺粒子の影響範囲計算可能粒子から圧力による反発力を受けて加速度を修正している
                                    double dist = sqrt(dist2);
                                    double w = WEI(dist, r);
                                    w *= (Prs[j] - pre_min) / dist2;
                                    Acc_x += v0 * w;
                                    Acc_y += v1 * w;
                                    Acc_z += v2 * w;
                                }
                            }
                            j = nxt[j];
                            if (j == -1) break;
                        }
                    }
                }
            }
            Acc[i * 3] = Acc_x * invDns[FLD] * A3;
            Acc[i * 3 + 1] = Acc_y * invDns[FLD] * A3;
            Acc[i * 3 + 2] = Acc_z * invDns[FLD] * A3;
        }
    }
}

//関数12 加速度に沿った移動二回目
void UpPcl2(void) {
    for (int i = 0; i < nP; i++) {
        if (Typ[i] == FLD) {
            Vel[i * 3] += Acc[i * 3] * DT;
            Vel[i * 3 + 1] += Acc[i * 3 + 1] * DT;
            Vel[i * 3 + 2] += Acc[i * 3 + 2] * DT;
            Pos[i * 3] += Acc[i * 3] * DT * DT;
            Pos[i * 3 + 1] += Acc[i * 3 + 1] * DT * DT;
            Pos[i * 3 + 2] += Acc[i * 3 + 2] * DT * DT;
            Acc[i * 3] = Acc[i * 3 + 1] = Acc[i * 3 + 2] = 0.0;
            ChkPcl(i);
        }
    }
}

void ClcEMPS(void) {
    while (1) {
        if (iLP % 100 == 0) {
            //標準出力
            int p_num = 0;
            for (int i = 0; i < nP; i++) { if (Typ[i] != GST)p_num++; }
            printf("%5d th TIM: %lf / p_num: %d\n", iLP, TIM, p_num);
        }
        if (iLP % OPT_FQC == 0) {
            //ファイル出力
            WrtDat();
            if (TIM >= FIN_TIM) { break; }
        }
        MkBkt();//バケット設定
        VscTrm();//粘性項重力項
        UpPcl1();//移動
        ChkCol();//剛体判定check collision
        MkPrs();//仮圧力
        PrsGrdTrm();//圧力勾配加速度
        UpPcl2();//移動
        MkPrs();//仮圧力
        /// TODO なんでMkPrsを二回呼ぶと二回目は陰解法として扱えるのか理解する
        for (int i = 0; i < nP; i++) { pav[i] += Prs[i]; }
        iLP++;
        TIM += DT;
    }
}

#include <sys/time.h>

double get_dtime(void) {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return ((double) (tv.tv_sec) + (double) (tv.tv_usec) * 0.000001);
}

int output(void);
int main(int argc, char **argv) {
    output();
    printf("start emps.\n");
    RdDat();
    AlcBkt();
    SetPara();

    double timer_sta = get_dtime();

    ClcEMPS();

    double timer_end = get_dtime();
    printf("Total        : %13.6lf sec\n", timer_end - timer_sta);


    free(Acc);
    free(Pos);
    free(Vel);
    free(Prs);
    free(pav);
    free(Typ);
    free(bfst);
    free(blst);
    free(nxt);
    printf("end emps.\n");
    return 0;
}
