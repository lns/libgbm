#ifndef DEFS_H
#define DEFS_H
/* 
 * Defines NAN,INF
 * INFINITY,POS_INF,NEG_INF
 */
#ifndef NAN
	#ifdef __INTEL_COMPILER
		#define NAN ((__builtin_nanl("")))
	#else
		#include <cmath>
		#ifndef NAN
		#define NAN (-sqrt(-1))
		#endif
	#endif
#endif

#ifndef INF
	#ifdef __INTEL_COMPILER
		#define INF ((__builtin_huge_vall()))
	#else
		#include <cmath>
		#ifdef INFINITY
			#define INF INFINITY
		#else
			#define INF ((double)1/0)
		#endif
	#endif
#endif

#ifndef INFINITY
#define INFINITY INF
#endif
#ifndef POS_INF
#define POS_INF INF
#endif
#ifndef NEG_INF
#define NEG_INF (-INF)
#endif
#endif

