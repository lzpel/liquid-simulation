#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "defs.h"
#include "bucket.h"

#define IN_FILE "./initCollapse.prof"
#define PCL_DST 0.02                    //平均粒子間距離
#define MIN_X  (0.0 - PCL_DST*3)    //解析領域のx方向の最小値
#define MIN_Y  (0.0 - PCL_DST*3)    //解析領域のy方向の最小値
#define MIN_Z  (0.0 - PCL_DST*3)    //解析領域のz方向の最小値
#define MAX_X  (1.0 + PCL_DST*3)    //解析領域のx方向の最大値
#define MAX_Y  (0.2 + PCL_DST*3)    //解析領域のy方向の最大値
#define MAX_Z  (0.6 + PCL_DST*30)    //解析領域のz方向の最大値

#define GST int(Type::GHOST)            //計算対象外粒子の種類番号
#define FLD int(Type::FLUID)                //流体粒子の種類番号
#define WLL int(Type::WALL)            //壁粒子の種類番号
#define NUM_TYP  10        //粒子の種類数
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
double r, r2;
Bucket bucket;
double n0, lmd, A1, A2, A3, rlim, rlim2, COL;
double Dns[NUM_TYP], invDns[NUM_TYP];

struct Particle {
    double Acc[3], Pos[3], Vel[3];
    double Prs, pav;
    int Typ;
};

std::vector<Particle> ps;

//関数01 計算領域を飛び出たら幽霊粒子に設定し計算から排除
void ChkPcl(Particle &p) {
    if (p.Pos[0] > MAX_X || p.Pos[0] < MIN_X ||
        p.Pos[1] > MAX_Y || p.Pos[1] < MIN_Y ||
        p.Pos[2] > MAX_Z || p.Pos[2] < MIN_Z) {
        p.Typ = GST;
        p.Prs = p.Vel[0] = p.Vel[1] = p.Vel[2] = 0.0;
    }
}

//関数02 初期状態読み込み
void RdDat(void) {
    fp = fopen(IN_FILE, "r");
    double distance, density;
    fscanf(fp, "%lf %lf", &distance, &density);
    ps = std::vector<Particle>();
    while (true) {
        int idx;
        Particle p = {};
        //粒子番号　粒子種別(GST -1 FLD 0 WLL) 座標三次元　速度三次元　（加速度は0）
        if (fscanf(fp, " %d %lf %lf %lf", &p.Typ, &p.Pos[0], &p.Pos[1], &p.Pos[2]) == EOF) {
            nP = ps.size();
            break;
        }
        if ((DIM == 2) && (p.Pos[1] != 0)) {
            continue;
        }
        ps.push_back(p);
    }
    printf("nP: %d\n", nP);
    fclose(fp);
    for (int i = 0; i < nP; i++) { ChkPcl(ps[i]); }
    for (int i = 0; i < nP; i++) { ps[i].Acc[0] = ps[i].Acc[1] = ps[i].Acc[2] = 0.0; }
}

//関数03 状態書き込み
void WrtDat(void) {
    char outout_filename[256];
    sprintf(outout_filename, "outputNew%05d.prof", iF);
    fp = fopen(outout_filename, "w");
    fprintf(fp, "%d\n", nP);
    for (int i = 0; i < nP; i++) {
        Particle &p = ps[i];
        fprintf(fp, "%d %lf %lf %lf %lf %lf %lf %lf %lf\n", p.Typ, p.Pos[0], p.Pos[1], p.Pos[2], p.Vel[0], p.Vel[1], p.Vel[2], p.Prs,
                p.pav / OPT_FQC);
        //書き込みのタイミングで時間平均圧力をリセットしている
        p.pav = 0.0;
    }
    fclose(fp);
    iF++;
}

//関数04 バケット構造のメモリを確保する
void AlcBkt(void) {
    //なんで2.1なの？
    r = PCL_DST * 2.1;        //影響半径
    r2 = r * r;
    bucket.min[0] = MIN_X;
    bucket.min[1] = MIN_Y;
    bucket.min[2] = MIN_Z;
    bucket.max[0] = MAX_X;
    bucket.max[1] = MAX_Y;
    bucket.max[2] = MAX_Z;
    bucket.rangeUpdateFinish(ps.size(), r * (1.0 + CRT_NUM));
}

//関数05 粒子法にかかわるパラメータの設定
void SetPara(void) {
    //理想的な格子の状態で粒子数密度とラプラシアンモデルの係数を求めている
    double tn0 = 0.0;
    double tlmd = 0.0;
    for (int ix = -4; ix < 5; ix++) {
        for (int iy = -4; iy < 5; iy++) {
            if (DIM == 2 && iy != 0)continue;
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
    n0 = tn0;            //p30 初期粒子数密度
    lmd = tlmd / tn0;    //p30 ラプラシアンモデルの係数λ

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
    bucket.clear();
    for (int i = 0; i < nP; i++) {
        Particle &p = ps[i];
        //計算しない
        if (p.Typ == GST)continue;
        bucket.add(p.Pos, i);
    }
}

//関数07 粘性項と重力を計算する関数
void VscTrm() {
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        if (pi.Typ == FLD) {
            double Acc_x = 0.0;
            double Acc_y = 0.0;
            double Acc_z = 0.0;
            double pos_ix = pi.Pos[0];
            double pos_iy = pi.Pos[1];
            double pos_iz = pi.Pos[2];
            double vec_ix = pi.Vel[0];
            double vec_iy = pi.Vel[1];
            double vec_iz = pi.Vel[2];
            for (int j = bucket.iterator(pi.Pos); j != -1; j = bucket.next(j)) {
                Particle &pj = ps[j];
                double v0 = pj.Pos[0] - pos_ix;
                double v1 = pj.Pos[1] - pos_iy;
                double v2 = pj.Pos[2] - pos_iz;
                double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                if (dist2 < r2) {
                    //近傍
                    if (j != i && pj.Typ != GST) {
                        //自分以外の流体粒子
                        //ウェイトと掛け算してAccに加算、Accの次元はまだ速度であるはず
                        double dist = sqrt(dist2);
                        double w = WEI(dist, r);
                        Acc_x += (pj.Vel[0] - vec_ix) * w;
                        Acc_y += (pj.Vel[1] - vec_iy) * w;
                        Acc_z += (pj.Vel[2] - vec_iz) * w;
                    }
                }
            }
            //係数A1によってAccの次元が速度から加速度に変わる、また重力を加算
            pi.Acc[0] = Acc_x * A1 + G_X;
            pi.Acc[1] = Acc_y * A1 + G_Y;
            pi.Acc[2] = Acc_z * A1 + G_Z;
        }
    }
}

//関数08 加速度の修正項から位置と速度を求める関数
void UpPcl1() {
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        if (pi.Typ == FLD) {
            pi.Vel[0] += pi.Acc[0] * DT;
            pi.Vel[1] += pi.Acc[1] * DT;
            pi.Vel[2] += pi.Acc[2] * DT;
            pi.Pos[0] += pi.Vel[0] * DT;
            pi.Pos[1] += pi.Vel[1] * DT;
            pi.Pos[2] += pi.Vel[2] * DT;
            pi.Acc[0] = pi.Acc[1] = pi.Acc[2] = 0.0;
            //計算領域を飛び出たら幽霊粒子に設定し計算から排除
            ChkPcl(pi);
        }
    }
}

//関数09 粒子の剛体衝突判定
void ChkCol() {
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        if (pi.Typ == FLD) {
            double mi = Dns[pi.Typ];
            double pos_ix = pi.Pos[0];
            double pos_iy = pi.Pos[1];
            double pos_iz = pi.Pos[2];
            double vec_ix = pi.Vel[0];
            double vec_iy = pi.Vel[1];
            double vec_iz = pi.Vel[2];
            double vec_ix2 = pi.Vel[0];
            double vec_iy2 = pi.Vel[1];
            double vec_iz2 = pi.Vel[2];
            for (int j = bucket.iterator(pi.Pos); j != -1; j = bucket.next(j)) {
                Particle &pj = ps[j];
                //相対座標v=xj-xi
                double v0 = pj.Pos[0] - pos_ix;
                double v1 = pj.Pos[1] - pos_iy;
                double v2 = pj.Pos[2] - pos_iz;
                double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                if (dist2 < rlim2) {
                    //衝突判定用二乗距離rlim2以内
                    if (j != i && pj.Typ != GST) {
                        //かつ自分以外の計算可能粒子
                        double fDT = (vec_ix - pj.Vel[0]) * v0 + (vec_iy - pj.Vel[1]) * v1 +
                                     (vec_iz - pj.Vel[2]) * v2;
                        if (fDT > 0.0) {
                            //veli-velj * xj-xi >0 で内積から対向を判定
                            //粒子種別に応じた密度
                            double mj = Dns[pj.Typ];
                            // TODO: この反発式は何？
                            fDT *= COL * mj / (mi + mj) / dist2;
                            vec_ix2 -= v0 * fDT;
                            vec_iy2 -= v1 * fDT;
                            vec_iz2 -= v2 * fDT;
                        }
                    }
                }
            }
            pi.Acc[0] = vec_ix2;
            pi.Acc[1] = vec_iy2;
            pi.Acc[2] = vec_iz2;
        }
    }
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        //TODO: 次元整合？
        pi.Vel[0] = pi.Acc[0];
        pi.Vel[1] = pi.Acc[1];
        pi.Vel[2] = pi.Acc[2];
    }
}

//関数10 仮の圧力を求める
void MkPrs() {
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        if (pi.Typ != GST) {
            double pos_ix = pi.Pos[0];
            double pos_iy = pi.Pos[1];
            double pos_iz = pi.Pos[2];
            double ni = 0.0;
            for (int j = bucket.iterator(pi.Pos); j != -1; j = bucket.next(j)) {
                Particle &pj = ps[j];
                //相対座標v=xj-xi
                double v0 = pj.Pos[0] - pos_ix;
                double v1 = pj.Pos[1] - pos_iy;
                double v2 = pj.Pos[2] - pos_iz;
                double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                if (dist2 < r2) {
                    if (j != i && pj.Typ != GST) {
                        double dist = sqrt(dist2);
                        double w = WEI(dist, r);
                        ni += w;
                    }
                }
            }
            double mi = Dns[pi.Typ];
            //技巧的な式、0か仮圧力が求まる、0は表面粒子、粒子毎にprs[i]に圧力を入れる
            double pressure = (ni > n0) * (ni - n0) * A2 * mi;
            pi.Prs = pressure;
        }
    }
}

//関数11 圧力勾配項を求める関数
void PrsGrdTrm() {
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        if (pi.Typ == FLD) {
            double Acc_x = 0.0;
            double Acc_y = 0.0;
            double Acc_z = 0.0;
            double pos_ix = pi.Pos[0];
            double pos_iy = pi.Pos[1];
            double pos_iz = pi.Pos[2];
            double pre_min = pi.Prs;
            for (int j = bucket.iterator(pi.Pos); j != -1; j = bucket.next(j)) {
                Particle &pj = ps[j];
                double v0 = pj.Pos[0] - pos_ix;
                double v1 = pj.Pos[1] - pos_iy;
                double v2 = pj.Pos[2] - pos_iz;
                double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                if (dist2 < r2) {
                    if (j != i && pj.Typ != GST) {
                        //自分以外の周辺粒子の影響範囲計算可能粒子に最低圧力を保証している
                        if (pre_min > pj.Prs)pre_min = pj.Prs;
                    }
                }
            }
            for (int j = bucket.iterator(pi.Pos); j != -1; j = bucket.next(j)) {
                Particle &pj = ps[j];
                double v0 = pj.Pos[0] - pos_ix;
                double v1 = pj.Pos[1] - pos_iy;
                double v2 = pj.Pos[2] - pos_iz;
                double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                if (dist2 < r2) {
                    if (j != i && pj.Typ != GST) {
                        //自分以外の周辺粒子の影響範囲計算可能粒子から圧力による反発力を受けて加速度を修正している
                        double dist = sqrt(dist2);
                        double w = WEI(dist, r);
                        w *= (pj.Prs - pre_min) / dist2;
                        Acc_x += v0 * w;
                        Acc_y += v1 * w;
                        Acc_z += v2 * w;
                    }
                }
            }
            pi.Acc[0] = Acc_x * invDns[FLD] * A3;
            pi.Acc[1] = Acc_y * invDns[FLD] * A3;
            pi.Acc[2] = Acc_z * invDns[FLD] * A3;
        }
    }
}

//関数12 加速度に沿った移動二回目
void UpPcl2(void) {
    for (int i = 0; i < nP; i++) {
        Particle &pi = ps[i];
        if (pi.Typ == FLD) {
            pi.Vel[0] += pi.Acc[0] * DT;
            pi.Vel[1] += pi.Acc[1] * DT;
            pi.Vel[2] += pi.Acc[2] * DT;
            pi.Pos[0] += pi.Acc[0] * DT * DT;
            pi.Pos[1] += pi.Acc[1] * DT * DT;
            pi.Pos[2] += pi.Acc[2] * DT * DT;
            pi.Acc[0] = pi.Acc[1] = pi.Acc[2] = 0.0;
            ChkPcl(pi);
        }
    }
}

void ClcEMPS(void) {
    while (1) {
        if (iLP % 100 == 0) {
            //標準出力
            int p_num = 0;
            for (int i = 0; i < nP; i++) { if (ps[i].Typ != GST)p_num++; }
            printf("%5d th TIM: %lf / p_num: %d\n", iLP, TIM, p_num);
            fflush(stdout);
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
        for (int i = 0; i < nP; i++) { ps[i].pav += ps[i].Prs; }
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

int main(int argc, char **argv) {
    printf("start emps.\n");
    RdDat();
    AlcBkt();
    SetPara();

    double timer_sta = get_dtime();

    ClcEMPS();

    double timer_end = get_dtime();
    printf("Total        : %13.6lf sec\n", timer_end - timer_sta);

    printf("end emps.\n");
    return 0;
}