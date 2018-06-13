#include <iostream>
#include <mpirxx.h>
#include "point.h"

int main() {
    mpz_t p, k;
    Point P, P1, P2, P3;
    mpz_init(k);

    mpz_init_set_str(p,"00FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFD97", 16);

    mpz_set_ui(P.x, 147ul);
    mpz_set_ui(P.z, 1ul);

    mpz_set_str(P1.x,"CBB8F5EBD80486B923EBFB17E5464173144CAC7B0447717B0EA8DE20545A6A23",16);
    mpz_set_ui(P1.z, 1ul);

    mpz_set_ui(P2.x, 345ul);
    mpz_set_ui(P2.z, 1ul);

    mpz_set_ui(P3.x, 1ul);
    mpz_set_ui(P3.z, 0ul);

    mpz_set_ui(k, 567ul);

    if (P.checkPoint(p) == 1)
        P.checkAnswer(p, k);
    else
        std::cout << "Точка " << P << " не лежит на эллиптической кривой.\n\n";

    if (P1.checkPoint(p) == 1)
        P1.checkAnswer(p, k);
    else
        std::cout << "Точка " << P1 << " не лежит на эллиптической кривой.\n\n";

    if (P2.checkPoint(p) == 1)
        P2.checkAnswer(p, k);
    else
        std::cout << "Точка " << P2 << " не лежит на эллиптической кривой.\n\n";

    P3.checkAnswer(p, k);

    mpz_clears(p, k, NULL);
    return 0;
}
