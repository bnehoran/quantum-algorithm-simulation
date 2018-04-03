/* 
 * main.c
 * --------
 * A simple main client to test unit functionality.
 * Only to be used during testing.
 * 
 * Author: Barak Nehoran
 * Advisor: Iasonas Petras
 */

#include "operator.h"
#include "state.h"
#include "vector.h"
#include "inverse.h"
#include <stdio.h>
#include <time.h>

int main(void) {
    srand(time(NULL));
    Vector_run_tests(100);
    printf("Vector tests: \t\tpassed.\n");
    State_run_tests(10);
    printf("State tests: \t\tpassed.\n");
    Operator_run_tests(10);
    printf("Operator tests: \tpassed.\n");
    Register_run_tests(100);
    printf("Register tests: \tpassed.\n");
    Inverse_run_tests();
    printf("All tests passed!\n");
    return(0);
}