//
// Created by quan on 4/13/17.
//

#include <utility.h>

double get_time()
{
    struct timeval t;
    double retval;

    gettimeofday(&t, NULL);
    retval = t.tv_sec;
    retval += t.tv_usec / 1000000.0;
    return retval;
}

