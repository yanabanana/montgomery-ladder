//
// Created by Yana Nemirovskaya on 13.12.17.
//

#ifndef MONTGOMERY_POINT_H
#define MONTGOMERY_POINT_H
#include <mpirxx.h>
#include <iostream>

struct Point {
    mpz_t x;
    mpz_t z;

    Point() {
        mpz_inits(x, z, NULL);
    }

    ~Point() {
        mpz_clears(x, z, NULL);
    }

    void plus(Point& P, const Point& P0, mpz_t p) {
        mpz_t m1, m2, m3, m4;
        mpz_inits(m1, m2, m3, m4, NULL);

        mpz_sub(m1, x, z);
        mpz_add(m2, P.x, P.z);
        mpz_mul(m1, m1, m2);

        mpz_add(m3, x, z);
        mpz_sub(m4, P.x, P.z);
        mpz_mul(m3, m3, m4);

        mpz_add(m2, m1, m3);
        mpz_pow_ui(m2, m2, 2ul);
        mpz_mul(m2, m2, P0.z);
        mpz_mod(x, m2, p);

        mpz_sub(m4, m1, m3);
        mpz_pow_ui(m4, m4, 2ul);
        mpz_mul(m4, m4, P0.x);
        mpz_mod(z, m4, p);

        mpz_clears(m1, m2, m3, m4, NULL);
    }


    void doubly(mpz_t p) {
        mpz_t m1, m2, m3, C;
        mpz_inits(m1, m2, m3, NULL);
        mpz_init_set_str(C,"2AF35AAD0B3DE2F60FC95AB8B9C2278C3CC58DC21B53907F53B3E54BCADB1BA1",16);

        mpz_sub(m1, x, z);
        mpz_pow_ui(m1, m1, 2ul);

        mpz_add(m2, x, z);
        mpz_pow_ui(m2, m2, 2ul);

        mpz_sub(m3, m2, m1);

        mpz_mul(m1, m1, m2);

        mpz_addmul(m2, C, m3);
        mpz_mul(m3, m3, m2);

        mpz_mod(x, m1, p);
        mpz_mod(z, m3, p);

        mpz_clears(m1, m2, m3, C, NULL);
    }

    void aff_coordinates(mpz_t p) {
        mpz_t invert;
        mpz_init(invert);
        if (mpz_cmp_ui(z, 0) != 0) {
            mpz_invert(invert, z, p);
            mpz_mul(x, x, invert);
            mpz_mod(x, x, p);
            mpz_mul(z, z, invert);
            mpz_mod(z, z, p);
        } else if (mpz_cmp_ui(x, 0) != 0) {
            mpz_invert(invert, x, p);
            mpz_mul(x, x, invert);
            mpz_mod(x, x, p);
        }
        mpz_clear(invert);
    }

    void power(const mpz_t k, mpz_t p) {
        Point R, P, R0, P0, neutral;
        mpz_set_ui(neutral.x, 1ul);
        mpz_set_ui(neutral.z, 0ul);
        if (*this == neutral)
            *this = neutral;
        else {
            int length;
            P0 = *this;
            R0 = P0;
            R = P0;
            R.doubly(p);
            length = mpz_sizeinbase(k, 2);
            for (int i = length - 2; i >= 0; i--) {
                if (!mpz_tstbit(k, i)) {
                    R.plus(*this, P0, p);
                    (*this).doubly(p);
                } else {
                    (*this).plus(R, R0, p);
                    R.doubly(p);
                }
            }
            (*this).aff_coordinates(p);
        }
    }

    int checkPoint(mpz_t p) {
        mpz_t A, B, y, a, x2, invB;
        int result;
        mpz_inits(a, y, x2, invB, NULL);
        mpz_init_set_str(A,"ABCD6AB42CF78BD83F256AE2E7089E30F31637086D4E41FD4ECF952F2B6C6E86",16);
        mpz_init_set_str(B,"ABCD6AB42CF78BD83F256AE2E7089E30F31637086D4E41FD4ECF952F2B6C6E88",16);

        mpz_pow_ui(a, x, 3ul);
        mpz_pow_ui(x2, x, 2ul);
        mpz_addmul(a, x2, A);
        mpz_add(a, a, x);
        mpz_invert(invB, B, p);
        mpz_mul(a, a, invB);
        mpz_mod(a, a, p);

        result = mpz_legendre(a, p);

        mpz_clears(A, B, y, a, x2, invB, NULL);
        return result;
    }


    void checkAnswer(mpz_t p, mpz_t k) {
        Point P1, P2, P3, P4, neutral;
        mpz_t q, q1, q2, m;
        mpz_init_set_str(q,"400000000000000000000000000000000FD8CDDFC87B6635C115AF556C360C67",16);
        mpz_init_set_ui(q1, 1ul);
        mpz_add(q1, q1, q);
        mpz_init_set_str(m, "01000000000000000000000000000000003F63377F21ED98D70456BD55B0D8319C", 16);
        P1 = *this;
        P2 = *this;
        P3 = *this;
        P4 = *this;
        mpz_set_ui(neutral.x, 1ul);
        mpz_set_ui(neutral.z, 0ul);
        P1.power(m, p);
        P2.power(q, p);
        P3.power(q1, p);
        P4.power(k, p);
        std::cout << "Исходная точка: " << *this << "\n";
        std::cout << "Исходная точка, возведенная в степень m - порядок группы точек эллиптической кривой: " << P1 << "\n";
        std::cout << "Исходная точка, возведенная в степень q - порядок циклической подгруппы точек эллиптической кривой: " << P2 << "\n";
        std::cout << "Исходная точка, возведенная в степень q + 1: " << P3 << "\n";
        std::cout << "Исходная точка, возведенная в степень k: " << P4 << "\n";
        if (P1 == neutral)
            std::cout << "Исходная точка в степени m равна нейтральному элементу => алгоритм возведения в степень работает корректно.\n";
        else
            std::cout << "Алгоритм возведения в степень работает некорректно.\n";
        if (P2 == neutral && P3 == *this)
            std::cout << "Исходная точка входит в циклическую подгруппу.\n\n\n";
        else
            std::cout << "Исходная точка не входит в циклическую подгруппу.\n\n\n";
        mpz_clears(q, q1, m, NULL);
    }


    Point& operator=(Point& P) {
        mpz_set(x, P.x);
        mpz_set(z, P.z);
        return *this;
    }


    bool operator==(Point& P) {
        if ((mpz_cmp(x, P.x) == 0) && (mpz_cmp(z, P.z) == 0))
            return true;
        else
            return false;
    }

    friend std::ostream& operator<<(std::ostream& output, Point& a) {
        output << "(" << a.x << "  ;  "<< a.z << ")";
        return output;
    }
};

#endif //MONTGOMERY_POINT_H


