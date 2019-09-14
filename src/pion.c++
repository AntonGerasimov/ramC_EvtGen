#include <iostream>
#include <fstream>
#include "ramC.h"
#include "./EvtWnPi2.hh"
#include "EvtGenBase/EvtPDL.hh"
using namespace std;


double mpi;

double spec_func_2pi(double q) {

    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi};
    double GF = 1;
    RAMBOC R(2, q, XM);
    double sum = 0;
    int nEv = 1e5;

    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 2 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        //       EvtVector4R q3 = R.getV(2);

        EvtVector4C alpha = cur.WCurrent(q1, q2);
        EvtVector4C beta = alpha.conj();
        EvtComplex gama = alpha * beta / (-3 * q * q);
        sum += real(gama) * wt;
        //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
        //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
            cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
            //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
            //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl;
        }
    };

    sum /= nEv;
    return sum;
}

double spec_func_3pi(double q) {

    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi, mpi};
    double GF = 1;
    RAMBOC R(3, q, XM);
    double sum = 0;
    int nEv = 1e5;

    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 3 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        EvtVector4R q3 = R.getV(2);

        EvtVector4C alpha = cur.WCurrent(q1, q2, q3);
        EvtVector4C beta = alpha.conj();
        EvtComplex gama = alpha * beta / (-3 * q * q);
        sum += real(gama) * wt;
        //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
        //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
            cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
            //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
            //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl;
        }
    };

    sum /= nEv;
    return sum;
}

double spec_func_4pi(double q) {

    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi, mpi, mpi};
    double GF = 1;
    RAMBOC R(4, q, XM);
    double sum = 0;
    int nEv = 1e5;

    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 4 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        EvtVector4R q3 = R.getV(2);
        EvtVector4R q4 = R.getV(3);

        EvtVector4C alpha = cur.WCurrent(q1, q2, q3, q4);
        EvtVector4C beta = alpha.conj();
        EvtComplex gama = alpha * beta / (-3 * q * q);
        sum += real(gama) * wt;
        //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
        //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
            cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
            //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
            //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl;
        }
    };

    sum /= nEv;
    return sum;
}

double spec_func_5pi(double q) {

    EvtWnPi2 cur;
    double pi = 3.14;
    double XM[] = {mpi, mpi, mpi, mpi, mpi};
    double GF = 1;
    RAMBOC R(5, q, XM);
    double sum = 0;
    int nEv = 1e5;

    for (int iEv = 0; iEv < nEv; ++iEv) {
        double wt = R.next();
        wt *= pow(2 * pi, 4 - 5 * 3);
        EvtVector4R q1 = R.getV(0);
        EvtVector4R q2 = R.getV(1);
        EvtVector4R q3 = R.getV(2);
        EvtVector4R q4 = R.getV(3);
        EvtVector4R q5 = R.getV(4);


        EvtVector4C alpha = cur.WCurrent(q1, q2, q3, q4, q5);
        EvtVector4C beta = alpha.conj();
        EvtComplex gama = alpha * beta / (-3 * q * q);
        sum += real(gama) * wt;
        //     double mtr2 = 128 * pow(GF, 2)*(p * q1)*(k * q1);
        //     sum += mtr2*wt;
        if (iEv < 0) {
            cout << " Debug print at iEv=" << iEv << " =========== " << endl;
            cout << " q1 = " << q1 << "; m2 = " << q1.mass2() << endl;
            cout << " q2 = " << q2 << "; m2 = " << q2.mass2() << endl;
            //           cout << " q3 = " << q3 << "; m2 = " << q3.mass2() << endl;
            //           cout << " q = " << q << "; m = " << p.mass() << endl;
            cout << "gama=" << gama << std::endl;
        }
    };

    sum /= nEv;
    return sum;
}

int save_sp(int n_pi) {
    string out_file_name = "./plot" + std::to_string(n_pi) + "pi.txt";
    cout << " Calculating spectral function for " << n_pi << " pions production. The results are saved to file " << out_file_name << endl;

    double qmax = n_pi;
    double qmin = n_pi*mpi;
    int qsize = 100;
    double qstep = (qmax - qmin) / qsize;


    ofstream out;
    out.open(out_file_name);

    for (int iq = 1; iq <= qsize; ++iq) {
        double qn = qmin + iq * qstep;
        double qsum;
        if (n_pi == 2)
            qsum = spec_func_2pi(qn);
        else if (n_pi == 3)
            qsum = spec_func_3pi(qn);
        else if (n_pi == 4)
            qsum = spec_func_4pi(qn);
        else if (n_pi == 5)
            qsum = spec_func_5pi(qn);
        else {
            cout << "n_pi = " << n_pi << " is not supported yet" << endl;
            return -1;
        };
        if (iq % (qsize / 10) == 0)
            cout << "========== " << 100. * iq / qsize << "% ===========" << endl;
        out << qn << " " << qsum << endl;
    };
    out.close();
    return 0;
}

int main(int argc, char *argv[]) {
    EvtPDL pdl;
    pdl.read("evt.pdl");
    mpi = EvtPDL::getMass(EvtPDL::getId("pi+"));

    if (argc != 2) {
        cout << " Wrong number of arguments! Use the format " << endl;
        cout << "\t  ./pion.exe <n_pi>" << endl;
        return -1;
    };
    int n_pi;
    if (string(argv[1]) == "all") {
        for (n_pi = 2; n_pi <= 5; ++n_pi) {
            save_sp(n_pi);
        };
        return 0;
    };

    n_pi = atoi(argv[1]);
    if (n_pi < 2 || n_pi > 5) {
        cout << " only 2 <= n_pi <=5 is supported" << endl;
        return -1;
    };

    save_sp(n_pi);


}

