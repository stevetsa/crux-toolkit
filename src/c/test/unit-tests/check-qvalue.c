#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "check-qvalue.h"
#include "q-value.h"
#include "crux-utils.h"

START_TEST (test_create){

  // A simple test of the decoy scoring routine.
  FLOAT_T scores1[] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};
  FLOAT_T scores2[] = {5.5, 4.5, 3.5, 2.5, 1.5, 0, -1, -2, -3, -4};
  FLOAT_T* qvalues = compute_decoy_qvalues(scores1, 10, scores2, 10);
  fail_unless(qvalues[0] == 0.0);
  fail_unless(qvalues[1] == 0.0);
  fail_unless(qvalues[2] == 0.0);
  fail_unless(qvalues[3] == 0.0);
  fail_unless(qvalues[4] == 0.0);
  fail_unless(qvalues[5] == (1.0 / 6.0));
  fail_unless(qvalues[6] == (2.0 / 7.0));
  fail_unless(qvalues[7] == (3.0 / 8.0));
  fail_unless(qvalues[8] == (4.0 / 9.0));
  fail_unless(qvalues[9] == (5.0 / 10.0));
  free(qvalues);

}
END_TEST

Suite *qvalue_suite(void){
  Suite *s = suite_create("qvalue_collection");
  TCase *tc_core = tcase_create("Core");
  suite_add_tcase(s, tc_core);
  tcase_add_test(tc_core, test_create);
  return s;
}
