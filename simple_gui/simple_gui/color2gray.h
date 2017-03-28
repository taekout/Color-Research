#ifndef COLOR2GRAY_H
#define COLOR2GRAY_H


#include <math.h>
#include <stdio.h>
#include <time.h>

#include "amy_colors.h"
#include "images.h"
#include "gui2.h"

#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_Shared_Image.H>
#include <FL/Fl.H>
#include <FL/Fl_PNM_Image.H>
#include <FL/Fl_Scroll.H>

 #define DEBUG_MODE 0

enum Imgi {I0,I45,I90,I135,I180,I225,I270,I315,IO,IOG,ING,INC};

extern Fl_PNM_Image *OutputImgs[12];


class Color2GrayIt{
 public :
  Color2GrayIt(int argc, char * argv[]);
  Color2GrayIt(char *name);
  ~Color2GrayIt();
  void init();
  int  parse_args(int argc, char * argv[]);
  void LoadSource();
  void run_c2g();
  void fillBoxWithImage(Fl_Box *box, int i, char *name);
  void DisplayCircle();
  void Display3();
  void Display4();
  void ShowSrcImg();
  void updateGUI();
  void getImageName(int i,  char *oname);

  void EnableQuantize(int v);

 void updateTheta(float t );

 void updateAlpha(float a );

 void updateMu(int m );
 void updateNumQColors(int n );
 void resizeSourceSm();
  
  //--------- Variables-------------
  float alpha;

  //variables controlling the quantization based acceleration.
  bool quantize;
  int q_colors;

  //theta is in radians.
  float theta;
  int r;

  char fname[256];
  char outname[256];
  char outname_color[256]; 
  char outname_OldGray[256];
  ColorImage source;  
  
  int SourceLoaded;
   float d2r;
  float r2d ;
Fl_Double_Window *win;
 Fl_Box *boxO, *boxNG, *boxNC, *boxG;
Fl_PNM_Image *OutputImgs[12];
 int WIN_CREATED;
 int NUM_OUT_IMGS;
};




#endif
