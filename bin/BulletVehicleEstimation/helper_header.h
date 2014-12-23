#ifndef HELPER_HEADER_H
#define HELPER_HEADER_H
//#include "heightfield1.h"//Replace whichever heightfield u want to include
#include "heightfield_random.h"
static char *MyHeightfield;
void inline creategreyscaledata()
{
  MyHeightfield = (char *)malloc(width*height);
  char pixel[3];//Pixel data
  char *pointer = header_data;
  int count;
#pragma omp parallel for private(count)
  for(count = 0;count < width*height; count++)
  {
    HEADER_PIXEL(pointer, pixel);
    MyHeightfield[count] = round(double(pixel[0] + pixel[1] + pixel[3])/3.0);
    //printf("Height: %d\n", MyHeightfield[count]);
  //  printf("Pixel: %d\t%d\t%d\n",-pixel[0],-pixel[1],-pixel[2]);
  }
  //delete header_data;
}
#endif
