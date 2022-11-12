#ifndef INCLUDED_UTILS_H
#define INCLUDED_UTILS_H

double degrees_to_radians(double degrees);

double random_double(void);
double random_double_in_range(double min, double max);
int random_int_in_range(int min, int max);

double clamp(double x, double min, double max);

#endif /* INCLUDED_UTILS_H */
