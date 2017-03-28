#ifndef IMAGES_H
#define IMAGES_H

#include <stdexcept>
#include <vector>
#include <algorithm>
#include "amy_colors.h"
  
using std::vector;
 


inline float sq(float s) { return s*s; }

struct ColorImage {

  amy_lab * data;
	
  typedef std::pair<amy_lab,int> amy_lab_int;
  vector<amy_lab_int> qdata;
	
  int colors;
	
  int w, h, N;

  ColorImage() : data(NULL) {}
  void clean() { if (data != NULL) delete [] data; }

  float calc_delta(int i, int j, float thetaL, float alphaL) const;
  float calc_qdelta(int i, int p, float thetaL, float alphaL) const;

  void load_quant_data(const char *fname);
	
  float * calc_d(float thetaL, float alphaL, int quantizeL );

  float * r_calc_d(int r, float thetaL, float alphaL);
	
  void load(const char * fname);
};

struct GrayImage {

  float * data;
  const int w, h, N;
  int i,j,k;

  //this will shift our data to best match the 
  //luminance channel of s.
  void post_solve(const ColorImage &s);

  GrayImage( ColorImage &s) : data(new float[s.N]), w(s.w), h(s.h), N(s.N) {
    for(i=0;i<N;i++) 
      data[i]=(s.data)[i].l;
  }
  ~GrayImage() { delete [] data; }

  void complete_solve(const float *d);
  void r_solve(const float *d, int r);
	
  void save(const char * fname) const;
  void saveColor(const char * fname, const ColorImage &source) const;
};	



#endif
