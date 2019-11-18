/* Program     : TOV Solver
 * File        : TOV_headers.h
 * Author      : Leo Werneck (werneck@if.usp.br)
 * Date        : October 29, 2019
 *
 * Description : This file contain the necessary function and variable
 *               definitions for the TOV Solver.
 *
 * Dependencies: None
 *
 * Reference(s): Read et al., PRD 79, 124032 (2009) | (https://arxiv.org/pdf/0812.2163.pdf)
 */
#ifndef __TOV_HEADERS_h__
#define __TOV_HEADERS_h__

#include <string>
using namespace std;

/* .---------.
 * | Defines |
 * .---------.
 */

/* Debugging option */
#define DEBUG 0

/* Definition of REAL */
#define REAL double

/* Physical constants */
#define SPEEDOFLIGHT 2.99792e+10
#define GNEWTON 6.67408e-08
#define SUNMASS 1.989e+33

/* Gridfunction definitions */
#define NGFS 3
#define PRESSURE 0
#define NU       1
#define MASS     2

/* EOS related definitions */
#define MAXEOSNAMESIZE 6
#define MAXEOSPARAMETERS 7

/* Integration related definitions */
#define MAXIMUMSTEPSIZE 1e-2
#define MAXINTEGRATIONS 300
#define MIN(a,b) a < b ? a : b

/* .------------.
 * | EOS Struct |
 * .------------.
 */
struct eos_struct
{
  int neos;
  REAL rho_b_PPEOS[MAXEOSPARAMETERS];
  REAL Gamma_PPEOS[MAXEOSPARAMETERS];
  REAL Kpoly_PPEOS[MAXEOSPARAMETERS];
  REAL Press_PPEOS[MAXEOSPARAMETERS];
  REAL epsIC_PPEOS[MAXEOSPARAMETERS];
};

#endif
