//
// Created by O'Connell, Jared on 7/26/17.
//

#ifndef AKT_MENDEL_H
#define AKT_MENDEL_H


/**
 * @name    duoMendel
 * @brief   returns true if duo genotypes are inconsistent with Mendel inheritance.
 *
 * @param [in] kid child  alternate allele count (3==missing)
 * @param [in] par parent alternate allele count (3==missing)
 *
 */
bool duoMendel(int kid, int par);

/**
 * @name    trioMendel
 * @brief   returns trio if duo genotypes are inconsistent with Mendel inheritance.
 *
 * @param [in] kid child  alternate allele count (3==missing)
 * @param [in] mum alternate allele count (3==missing)
 * @param [in] dad alternate allele count (3==missing)
 *
 */
bool trioMendel(int kid, int mum, int dad);

#endif //AKT_MENDEL_H
