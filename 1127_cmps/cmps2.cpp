//
// Created by misumi on 2021/11/27.
//

#include "time.h"
#include "stdio.h"
#include "math.h"
#include "vector"
#include "bucket.h"
#include "defs.h"

struct Particle {
    double Acc[3], Pos[3], Vel[3];
    double Prs, PrsSum;
    Type Typ;
    int nxt;//同バケット内の次の粒子
};

//定数
const int Dimention = 3;
const double InfluenceRadiusCoefficient = 2.1;
const double CourantConditionCoefficient = 0.1; //TODO: 意味を理解する
const double dt = 0.0005;//TODO: 時間刻み幅の安定条件を理解する
const double KNM_VSC = 0.000001;   //動粘性係数
const double dtWrite = 0.1;//出力時間刻み幅
const double G_X = 0.0;           //重力加速度のx成分
const double G_Y = 0.0;           //重力加速度のy成分
const double G_Z = -9.8;            //重力加速度のz成分
//値
double particleDistance;
double density;//kg/m^3 水なら1000
std::vector<Particle> ps;
Bucket bucket;
double influenceRadius;
double range[2][3];
double n0, lambda;

void readProf(const char *fname) {
    FILE *fp = fopen(fname, "r");
    fscanf(fp, "%lf %lf", &particleDistance, &density);
    ps = std::vector<Particle>();
    Particle p = {};
    while (fscanf(fp, " %d %lf %lf %lf", &p.Typ, &p.Pos[0], &p.Pos[1], &p.Pos[2]) != EOF) {
        if (Dimention == 2 && p.Pos[1] != 0) continue;
        bucket.rangeUpdate(p.Pos);
        ps.push_back(p);
    }
    fclose(fp);
    //バケット初期化
    influenceRadius = particleDistance * InfluenceRadiusCoefficient;
    printf("particleDistance = %lf, influenceRadius = %lf, ps.size() = %d\n", particleDistance, influenceRadius, ps.size());
    bucket.rangeUpdateFinishWithLength(influenceRadius * (1.0 + CourantConditionCoefficient));
    //粒子数密度
    double tn0 = 0.0, tlambda = 0.0;
    for (int ix = -4; ix < 5; ix++) {
        for (int iy = -4; iy < 5; iy++) {
            for (int iz = -4; iz < 5; iz++) {
                if (Dimention == 2 && iy != 0)continue;
                const double pos[3] = VectorScalar(ix, iy, iz, particleDistance);
                const double dist2 = InnerProduct(pos, pos);
                if (dist2 <= influenceRadius * influenceRadius) {
                    if (dist2 == 0.0)continue;
                    tn0 += Kernel(sqrt(dist2), influenceRadius);
                    tlambda += dist2 * Kernel(sqrt(dist2), influenceRadius);
                }
            }
        }
    }
    n0 = tn0;            //p30 初期粒子数密度
    lambda = tlambda / tn0;    //p30 ラプラシアンモデルの係数λ
    printf("n0 = %lf, lambda = %lf\n", n0, lambda);
}

void WriteProf() {
    static int iF = 0;
    char outout_filename[256];
    sprintf(outout_filename, __FILE__".%05d.prof", iF++);
    FILE *fp = fopen(outout_filename, "w");
    fprintf(fp, "%d\n", ps.size());
    for (int i = 0; i < ps.size(); i++) {
        Particle &p = ps[i];
        fprintf(fp, " %d %d %lf %lf %lf %lf %lf %lf %lf %lf\n", i, p.Typ, p.Pos[0], p.Pos[1], p.Pos[2], p.Vel[0], p.Vel[1], p.Vel[2], p.Prs, p.PrsSum * dt / dtWrite);
        p.PrsSum = 0.0;//書き込みのタイミングで時間平均圧力をリセットしている
    }
    fclose(fp);
}

void BucketUpdate() {
    bucket.clear();
    for (int i = 0; i < ps.size(); i++) { ps[i].nxt = -1; }
    for (int i = 0; i < ps.size(); i++) {
        Particle &p = ps[i];
        if (p.Typ == Type::GHOST)continue;
        p.nxt = bucket.getFirst(p.Pos);
        bucket.setFirst(p.Pos, i);
    }
}

void VscTrm() {
    for (int i = 0; i < ps.size(); i++) {
        Particle &pi = ps[i];
        if (pi.Typ == Type::FLUID) {
            double Acc_x = 0.0;
            double Acc_y = 0.0;
            double Acc_z = 0.0;
            double pos_ix = pi.Pos[0];
            double pos_iy = pi.Pos[1];
            double pos_iz = pi.Pos[2];
            double vec_ix = pi.Vel[0];
            double vec_iy = pi.Vel[1];
            double vec_iz = pi.Vel[2];
            const int *firsts = bucket.getFirsts(pi.Pos);
            for (int f = 0; f < 27; f++) {
                for (int j = firsts[f]; j != -1;) {
                    Particle &pj = ps[j];
                    //ここで周辺粒子が計算される
                    double v0 = pj.Pos[0] - pos_ix;
                    double v1 = pj.Pos[1] - pos_iy;
                    double v2 = pj.Pos[2] - pos_iz;
                    double dist2 = v0 * v0 + v1 * v1 + v2 * v2;
                    if (dist2 < influenceRadius * influenceRadius) {
                        //近傍
                        if (j != i && pj.Typ != Type::GHOST) {
                            //自分以外の流体粒子
                            //ウェイトと掛け算してAccに加算、Accの次元はまだ速度であるはず
                            double w = Kernel(sqrt(dist2), influenceRadius);
                            Acc_x += (pj.Vel[0] - vec_ix) * w;
                            Acc_y += (pj.Vel[1] - vec_iy) * w;
                            Acc_z += (pj.Vel[2] - vec_iz) * w;
                        }
                    }
                    j = pj.nxt;
                }
            }
            //係数A1によってAccの次元が速度から加速度に変わる、また重力を加算
            const double A1=2.0 * KNM_VSC * Dimention / n0 / lambda;
            pi.Acc[0] = Acc_x * A1 + G_X;
            pi.Acc[1] = Acc_y * A1 + G_Y;
            pi.Acc[2] = Acc_z * A1 + G_Z;
        }
    }
}

int main() {
    clockStart()
    readProf("./init20cm.prof");
    double time;
    for (int iter = 0; (time = iter * dt) < 1.0; iter++) {
        if (int(time / 0.1 + 1) != int(time / 0.1 + 1 - dt)) {
            //0.1秒毎に出力
            WriteProf();
            int p_num = 0;
            for (int i = 0; i < ps.size(); i++) { if (ps[i].Typ != Type::GHOST)p_num++; }
            printf("iter = %5d, time = %lf, validParticleNum = %d\n", iter, time, p_num);
        }
        BucketUpdate();//バケット設定
        VscTrm();//粘性項重力項
        //UpPcl1();//移動
        //ChkCol();//剛体判定check collision
        //MkPrs();//仮圧力
        //PrsGrdTrm();//圧力勾配加速度
        //UpPcl2();//移動
        //MkPrs();//仮圧力
        /// TODO なんでMkPrsを二回呼ぶと二回目は陰解法として扱えるのか理解する
        for (int i = 0; i < ps.size(); i++) { ps[i].PrsSum += ps[i].Prs; }
    }
    clockPrint()
}
