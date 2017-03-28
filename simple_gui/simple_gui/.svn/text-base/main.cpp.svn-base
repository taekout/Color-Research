///\file main.cpp
///minimalistic implementation of the color2gray algorithms.

//in an attempt to make this code as widely useable as possible, while minimizing
//dependencies, i'm using ppm images.

//Copyright 2005 sven olsen and amy gooch

#define _USE_MATH_DEFINES

#include "gui2.h"
#include "color2gray.h"
#include "thetacircle2.h"
#include "file_chooser.h"
#include <FL/Fl_Window.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Box.H>
#include <FL/Fl_Image.H>
#include <FL/Fl_Shared_Image.H>
#include <FL/Fl.H>
#include <FL/Fl_PNM_Image.H>

	
int win_w_size;
int win_h_size;

Color2GrayIt *c2g;
Color2GrayIt *c2gTC;

Fl_Double_Window *w_ctrl;

void EnableQuantize(Fl_Check_Button* b, void*){
  #if DEBUG_MODE
   printf("EnableQuantize button pressed with value  %d\n",
	 (int)b->value());
 #endif
    c2g->EnableQuantize((int)b->value());
}

void updateTheta(Fl_Dial* d, void*){
  #if DEBUG_MODE
  printf("updateTheta Dial pressed with value  %f\n",
	 (float)d->value());
 #endif
    c2g->updateTheta((float)d->value());
}
void updateThetaV(Fl_Value_Input* d, void*){
 #if DEBUG_MODE
   printf("updateTheta Input pressed with value  %f\n",
	 (float)d->value());
  #endif
   c2g->updateTheta((float)d->value());
}
 
void updateAlpha(Fl_Roller* r, void*){
 #if DEBUG_MODE
   printf("updateAlpha Roller pressed with value  %f\n",
	 (float)r->value());
  #endif
   c2g->updateAlpha((float)r->value());
}
void updateAlphaV(Fl_Value_Input* r, void*){
  #if DEBUG_MODE
  printf("updateAlphaV Input pressed with value  %f\n",
	 (float)r->value());
  #endif
   c2g->updateAlpha((float)r->value());
}
 

void updateMu(Fl_Roller* r, void*){
  #if DEBUG_MODE
  printf("updateMu Roller pressed with value  %d\n",
	(int)r->value());
 #endif
    c2g->updateMu((int)r->value());
}
void updateMuV(Fl_Value_Input* r, void*){
  #if DEBUG_MODE
  printf("updateMu Input pressed with value  %d\n",
	 (int)r->value());
  #endif
   c2g->updateMu((int)r->value());
}
 
void updateNumQColors(Fl_Value_Input* i, void*){
 #if DEBUG_MODE
   printf("NumQColors Input pressed with value  %d\n",
	 (int)i->value());
  #endif
   c2g->updateNumQColors((int)i->value());
}
void LoadSource(Fl_Button *b, void *){
  CreateFileBrowser();
}

void run_c2g_thetaCircle(Fl_Button *b, void *){
  strcpy(c2gTC->fname, c2g->fname);
  c2gTC->resizeSourceSm();  //Make source image fit 200
  c2gTC->q_colors = 30;
  c2gTC->quantize = true;
  c2gTC->r = 0;
  for(int t = 0; t<360; t=t+45){
    c2gTC->theta = t*c2gTC->d2r;
    c2gTC->run_c2g();
  }
  
  c2gTC->DisplayCircle();
}

void run_c2g_callback(Fl_Button *b, void *){
  c2g->run_c2g();
  //c2g->Display4();
  c2g->Display4();
}

// another callback to quit the program
void quit_cb(Fl_Button* , void *)
{
  exit(1);
}

void ShowSrcImg(Fl_Box*, void*){
  c2g->ShowSrcImg();
 
}
 

int main(int argc, char * argv[]) {

  c2g = new Color2GrayIt(argc,argv);

  if (argc >=2)
    c2gTC = new Color2GrayIt(argv[1]);
  else 
    c2gTC = new Color2GrayIt("cb45sm.ppm");

  fl_register_images();

  w_ctrl =  make_ctrl_window();
  //  Fl_Double_Window *w_ctrl = thetaCircleWindow();
  (*w_ctrl).resizable();
  (*w_ctrl).end();
  (*w_ctrl).show();
  
  
  c2g->updateGUI();
  printf("redrawing ctrl window\n");
  (*w_ctrl).redraw();
  
  //This has to be the last thing in this function
  return(  Fl::run());
}

