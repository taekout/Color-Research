
#include "color2gray.h"
#include "thetacircle2.h"

extern Color2GrayIt *c2g;

inline float crunch(float chrom_dist, float alphaL) {
  return alphaL == 0 ? 0 :  alphaL*tanh(chrom_dist/alphaL);
}

float * ColorImage::calc_d(float thetaL, float alphaL, int quantizeL) {
  float * d = new float [N];
  int i,j;
  for( i=0; i<N;i++) d[i] = 0;

  if(quantizeL) {
    printf("theta = %f (%f deg) alpha = %f\n", thetaL, thetaL*180/M_PI,alphaL);
    int p;
    for( i=0; i<N;i++) for(p=0;p<qdata.size();p++)
      d[i]+=calc_qdelta(i,p,thetaL,alphaL);
		
  }
  else {

    //more obvious but slower code for the unquantized full solve.
    /*for( i=0; i<N;i++) for(j=0;j<N;j++) {
      float delta = calc_delta(i,j);
      d[i]+=delta;
      }*/

    for( i=0; i<N;i++) for(j=i+1;j<N;j++) {
      float delta = calc_delta(i,j,thetaL,alphaL);
      d[i]+=delta;
      d[j]-=delta;
    }
    
  }
  return d;
}

float * ColorImage::r_calc_d(int r, float thetaL, float alphaL) {
  float * d = new float [N];
  int i;
  for( i=0; i<N;i++) d[i] = 0;

  int x, y;
  for(x=0;x<w;x++) for(y=0;y<h;y++) {
    int xx,yy;

    i=x+y*w;

    for(xx=x-r;xx<=x+r;xx++) {
      if(xx<0 || xx>=w) continue;
      for(yy=y-r;yy<=y+r;yy++) 		{
	if (yy>=h || yy<0) continue;	
	int j=xx+yy*w;
	float delta = calc_delta(i,j, thetaL, alphaL);
	d[i]+=delta;
	d[j]-=delta;				
      }
    }
  }
  return d;
}

void GrayImage::r_solve(const float *d, int r) {

  const int iters=30;
  int k, x, y;

  for(k=0;k<iters;k++) {

    printf("iter %d\n.",k);

    //perform a Gauss-Seidel relaxation.
    for(x=0;x<w;x++) for(y=0;y<h;y++) {
      float sum=0;
      int count=0;
      int xx,yy;

      for(xx=x-r;xx<=x+r;xx++) {
	if(xx<0 || xx>=w) continue;
	for(yy=y-r;yy<=y+r;yy++) 		{
	  if (yy>=h || yy<0) continue;	
	  sum+=data[xx+yy*w];
	  count++;
	}
      }
			
      data[x+y*w]=(d[x+w*y] +  sum) / (float)count;
      
    }
  }
}



void GrayImage::complete_solve(const float * d) {
  for(i=1;i<N;i++) {
    data[i] = d[i]-d[i-1]+N*data[i-1];
    data[i] /= (float)N;
  }
}

void GrayImage::post_solve(const ColorImage &s) {
  float error=0;
  for(i=0;i<N;i++) 
    error+=data[i]-(s.data)[i].l;
  error/=N;
  for( i=0;i<N;i++) data[i]=data[i]-error;
}


bool qchoice(char c, int *quantizeL, int *q_colors) {
  if(c=='y') {
    (*quantizeL)=true;
    printf("\nuse how many colors (64 is usually a sensible choice)?\n");
    scanf("%d",(q_colors));
  }
  else if(c=='n') {
    (*quantizeL)=false;
  }
  else {
    printf("\nplease enter 'y' or 'n'.\n");
    return false;
  }
  return true;
}


float ColorImage::calc_delta(int i, int j, float thetaL, float alphaL) const {
  const amy_lab &a = data[i];
  const amy_lab &b = data[j];

  float dL=a.l-b.l;
  float dC= crunch(sqrt(sq(a.a-b.a)+sq(a.b-b.b)), alphaL);

  if(fabsf(dL)>dC) return dL;
  return dC*(  (cos(thetaL)*(a.a-b.a) + sin(thetaL)*(a.b-b.b)) > 0 ? 1 : -1);
}


float ColorImage::calc_qdelta(int i, int p, float thetaL, float alphaL) const {
	const amy_lab &a = data[i];
	const amy_lab &b = qdata[p].first;

	float dL=a.l-b.l;
	float dC= crunch(sqrt(sq(a.a-b.a)+sq(a.b-b.b)),alphaL);

	if(fabsf(dL)>dC) return qdata[p].second*dL;
	return qdata[p].second*dC*(  (cos(thetaL)*(a.a-b.a) + sin(thetaL)*(a.b-b.b)) > 0 ? 1 : -1);
}
 

Color2GrayIt::~Color2GrayIt(){
  source.clean();
}

// aaa
void Color2GrayIt::init(){
 alpha = 10;
  quantize = true;
  q_colors = 64;
  theta = 3.14/4.0;
  r = 0;
  d2r = 3.14159 / 180.0;
  r2d = 180.0 / 3.14159;
  SourceLoaded = 0;
  WIN_CREATED = 0;
  win = NULL;
   int NUM_OUT_IMGS = 12;
   //   *OutputImgs = (Fl_PNM_Image *)malloc(size(Fl_PNM_Image *)*NUM_OUT_IMGS);

  for(int i=0; i<NUM_OUT_IMGS; i++)
    OutputImgs[i] = NULL;
}
Color2GrayIt::Color2GrayIt(int argc, char * argv[]){
  init();
  parse_args(argc,argv);
 
}

Color2GrayIt::Color2GrayIt(char *name){
  printf("init color2gray: %s\n", name);
  init();  strcpy(fname, name);

  
}

void fillBoxWithImage(Fl_Box *box, int i, char *name, Color2GrayIt *Lc2g){
  Lc2g->fillBoxWithImage(box,i,name);
}


void Color2GrayIt::fillBoxWithImage(Fl_Box *box, int i, char *name){
  if (OutputImgs[i] != NULL){
    delete(OutputImgs[i]);
  }
  OutputImgs[i] = new Fl_PNM_Image(name); 
  (*box).image(*(OutputImgs[i]));
}

#define MAX(a,b) a < b ? b : a

void Color2GrayIt::LoadSource(){

  SourceLoaded = 1;
  printf("Loading source image %s for colo2gray\n", fname);
  //load the image
  //ColorImage source;
  source.load(fname);

  
  //GrayImage dest(source);	
//   Fl_Window w(100,100,source.w,source.h, "Source Image");
//   (w).begin();
 
//   Fl_PNM_Image Oimg(fname);

//   // Fl_Shared_Image sfl_img(fname, flimg);
//   Fl_Box boxO(0,0,source.w,source.h);

//   boxO.image(Oimg);

//   (w).end();
//   (w).show();
//   //Fl::run();
}

int Color2GrayIt::parse_args(int argc, char **argv){

  if (argc >=2)
    strcpy(fname,argv[1]);
  else
    strcpy(fname,"cb45sm.ppm");

  for(int i=2; i< argc; i++){
    if (!strcmp(argv[i], "-theta")){
      i++;
      sscanf(argv[i], "%f", (&theta));
      theta=theta*d2r ;
      printf("Theta = %.2f\n", theta*180/M_PI);
    }
    else if (!strcmp(argv[i], "-alpha")){;
      i++;
      sscanf(argv[i], "%f", &alpha );
      printf("Alpha = %.2f\n", alpha);
    }
    else if (!strcmp(argv[i], "-r")){
      i++;
      sscanf(argv[i], "%d", &r );
      printf("mu = %d\n", r);
    }
    else if (!strcmp(argv[i], "-mu")){
      i++;
      sscanf(argv[i], "%d", &r );
      printf("mu = %d\n", r);
    }
    else  if (!strcmp(argv[i], "-q")){
      i++;
      quantize = true;
      sscanf(argv[i], "%d", &q_colors );
      printf("q = %d\n", q_colors);
    }
    else{
      i++; printf("%s is not a valid option: color2gray file.ppm -theta 45 -alpha 10 -r 0 -q 0\n", argv[i]);
    }

  }//end of for i on parse args
  fprintf(stderr, "Done with arg parsing\n");




  //print out an informative message.
  //  printf("Arguments: color2gray algorithm on \"%s\"\n with alpha=%f, theta=%f ", 
  //	 fname, alpha,theta*r2d);
  if(r==0) {
    printf("(complete case.)\n");
    //     if (!PARSE_ARGS_AMY){
    //       printf("\nwould you like to use the quantization based inexact version of the algorithm?\n (y/n)? \t(requires image magick's 'convert' in the path.)\n");
    //       char c; scanf("%c",&c);
    //       while(!qchoice(c)) scanf("%c",&c);
    //     }
  }
  else printf(", r=%d\n",r);



  return(1);
}

//Use Image magick to resize image to fit in 200x200 for theta circle & testing
void Color2GrayIt::resizeSourceSm(){
  int newSize = 200;
  char m_string[256];
  if ( SourceLoaded ) {
    source.clean();
  }

  printf("resizeSourceSm: loading image.. \n");  
  LoadSource();  
  
  if (source.w > 200 || source.h > 200){
    sprintf(m_string,"convert -resize %dx%d %s %d_%s",newSize,newSize,fname, newSize,fname);
    printf("\nUsing image magick to create a reduced image size: \n%s\n",m_string);
    system(m_string);
    printf("done.\n\n");
    sprintf(m_string,"%d_%s",newSize,fname);
    strcpy(fname, m_string);
    //remove old & reset
    source.clean();
    SourceLoaded = 0;
  }
}


void Color2GrayIt::updateGUI(){
  thetaD->value(theta*r2d);
  alphaR->value(alpha);
  muR->value(r);
  alphaV->value(alpha);
  muV->value( r);
  thetaV->value (theta*r2d);
  numQColorsV->value(q_colors);
  
  if (win != NULL)
    (*win).redraw();
  extern Fl_Double_Window *w_ctrl;
  (*w_ctrl).redraw();
}

void Color2GrayIt::run_c2g(){
  printf("Running Color2Gray....\n");
  if ( !SourceLoaded ){
    printf("Will call LoadSource with image %s\n", fname);
    LoadSource();  

  }
  else     printf("Assuming %s is already loaded\n", fname);
  // Create a meaningful output file name for the color2gray image:
  // Note: this only works if the filename is just the file and doesn't include a directory path
  if (r == 0)
    sprintf(outname, "c2g_theta%.1f_a%.1f_mu%s_q%d.%d_%s",theta*r2d,alpha,"FULL",quantize,q_colors,fname);
  else
    sprintf(outname, "c2g_theta%.1f_a%.1f_mu%d_q%d.%d_%s",theta*r2d,alpha,r,quantize,q_colors,fname);

  if (r == 0)
    sprintf(outname_color, "c2gC_theta%.1f_a%.1f_mu%s_q%d.%d_%s",theta*r2d,alpha,"FULL",quantize,q_colors,fname);
  else
    sprintf(outname_color, "c2gC_theta%.1f_a%.1f_mu%d_q%d.%d_%s",theta*r2d,alpha,r,quantize,q_colors,fname);
  printf("Color2Gray result will be output to: \n\t%s\n and \t%s\n", outname, outname_color); 

 

   printf("Running Color2Gray with:\n\tTheta=%f deg.\n\tAlpha=%f\n\tMu=%d\n\tQuantize=%d\n\tNumColors=%d\n\t",
	  theta*r2d,alpha,r,quantize,q_colors);

  clock_t start = clock();
  GrayImage dest(source);

  //output original gray image:
  sprintf(outname_OldGray, "OldGray_%s",fname);
  dest.save(outname_OldGray);

  //solve, either using the complete case or the neighboorhod case.
  float *d;
  if(r) { 
    d = source.r_calc_d(r, theta, alpha);
    dest.r_solve(d,r);
  }
  else {
    if(quantize) {
      char magick_string[256];
      sprintf(magick_string,"convert -colors %d %s quant-%s",q_colors,fname,fname);
      printf("\nUsing image magick to create a quantized image: \n%s\n",magick_string);
      system(magick_string);
      printf("done.\n\n");
      sprintf(magick_string,"quant-%s",fname);
      source.load_quant_data(magick_string);
    }
	
	
    d = source.calc_d(theta,alpha, quantize);
    dest.complete_solve(d);
  }
  dest.post_solve(source);

  clock_t end = clock();

  printf("\nc2g completed in %.03f seconds.\n",(end-start)/(float)CLOCKS_PER_SEC);

  dest.save(outname);
  dest.saveColor(outname_color, source);

  delete [] d;
  //  source.clean();
  if (WIN_CREATED > 0){
    fillBoxWithImage(boxNG,ING,outname);
    fillBoxWithImage(boxNC,INC,outname_color);
    win->redraw();
  }
  
}

void drawImg(){
  
}

void Color2GrayIt::getImageName(int i, char *oname){
  
  printf("Original file name = %s\n", fname);
  switch (i){
  case -2://gray image
    printf("getImageName %d -2 case\n", i);
    strcpy(oname, outname_OldGray);
    break;
  case -1://gray image
    printf("getImageName %d -1 case\n",i);
    strcpy(oname, fname);
    break;
  default://gray image
    printf("getImageName %d default case\n",i);
    if (r == 0)
      sprintf(oname, "c2g_theta%.1f_a%.1f_mu%s_q%d.%d_%s",(float)i,alpha,"FULL",quantize,q_colors,fname);
    else
      sprintf(oname, "c2g_theta%.1f_a%.1f_mu%d_q%d.%d_%s",(float)i,alpha,r,quantize,q_colors,fname);
   
    break;
   
  }
  printf("getImageName %d = %s\n",i, oname);
  
}
void Color2GrayIt::DisplayCircle(){
  fprintf(stderr, "Should be updating NewGray and New Color images over various theta values....\n");
 
  win = thetaCircleWindow();
  
  //  (*wTC).end();
  (*win).show();
  Fl::run();
}



void Color2GrayIt::Display3(){
  if (WIN_CREATED==0) {
    WIN_CREATED = 3;
    fprintf(stderr, "Should be updating NewGray and New Color images....\n");
    fl_register_images();
    win = new Fl_Double_Window( source.w*3+100, source.h*1+100, "Color2Gray Results");
    (*win).begin();
    (*win).color((Fl_Color)55);
    Fl_Scroll *scr; 
    if ((source.w*3+100) > 1440 || (source.h*1+100) > 900){
      printf("Using scrollscreen: window would be %d x %d\n", source.w*3+100, source.h*1+100);
      (*win).resize(0,0,810,810);
    scr= new Fl_Scroll(0,0,800,800);
    Fl_Group::current()->resizable(scr);
    }
    Fl_Group::current()->resizable(*win);
    //   // Fl_Shared_Image sfl_img(fname, flimg);
    boxO = new Fl_Box(50,50,source.w,source.h, "Original");
    boxNG =  new Fl_Box(source.w+10+50,50,source.w,source.h, "New Gray");
    boxNC =  new Fl_Box((source.w*2)+50+10+10,50,source.w,source.h, "New Color");
    fillBoxWithImage(boxO,IO,fname);
    fillBoxWithImage(boxNG,ING,outname);
    fillBoxWithImage(boxNC,INC,outname_color);
  
    (*win).resizable();
    (*win).end();
    
    (*win).show();
    Fl::run();  
  }
  else {
    printf("Redraw Display3()\n");
    fillBoxWithImage(boxNG,ING,outname);
    fillBoxWithImage(boxNC,INC,outname_color); 
    win->redraw();
 }
}

void Color2GrayIt::Display4(){
  if (WIN_CREATED==0){
    WIN_CREATED = 4;
    fprintf(stderr, "Should be updating NewGray and New Color images....\n");
    fl_register_images();
    win = new Fl_Double_Window( source.w*2+100, source.h*2+100, "Color2Gray Results");
    (*win).begin();
    (*win).color((Fl_Color)55);
    Fl_Scroll *scr; 
    if ((source.w*2+100) > 1440 || (source.h*2+100) > 900){
    printf("Using scrollscreen: window would be %d x %d\n", source.w*2+100, source.h*2+100);
    (*win).resize(0,0,810,810);
    scr= new Fl_Scroll(0,0,800,800);
    Fl_Group::current()->resizable(scr);
    }
    Fl_Group::current()->resizable(*win);		\
  //   // Fl_Shared_Image sfl_img(fname, flimg);
    boxO =  new Fl_Box(50,30,source.w,source.h, "Original");
    boxG = new Fl_Box(source.w+50+10,30,source.w,source.h, "Old Gray");
    boxNG  = new Fl_Box(50,30+40+source.h,source.w,source.h, "New Gray");
    boxNC  =  new Fl_Box((source.w)+50+10,30+40+source.h,source.w,source.h, "New Color");
    fillBoxWithImage(boxO,IO,fname);
    fillBoxWithImage(boxG,IOG,outname_OldGray);
    fillBoxWithImage(boxNG,ING,outname);
    fillBoxWithImage(boxNC,INC,outname_color);
    
    (*win).resizable();
    (*win).end();
    
    (*win).show();
    Fl::run();  
  }
  else {
    printf("Redraw Display4()\n");
    fillBoxWithImage(boxNG,ING,outname);
    fillBoxWithImage(boxNC,INC,outname_color);
    win->redraw();
  }
}
void Color2GrayIt::ShowSrcImg(){
  Fl_PNM_Image Oimg(fname);
  
  //Fl_Shared_Image sfl_img(fname, flimg);
  boxO = new Fl_Box(0,0,source.w,source.h);  
  (*boxO).image(Oimg);
}


void Color2GrayIt::EnableQuantize(int v){
#if DEBUG_MODE
  printf("EnableQuantize Function called with value %d\n", v);
#endif
  quantize = (bool)v;
  updateGUI();
}

void Color2GrayIt::updateTheta(float t){
#if DEBUG_MODE
  printf("updateTheta Function called with value %f\n", t);
#endif
  theta = t*d2r;
  updateGUI();
  
}

void Color2GrayIt::updateAlpha(float a ){
#if DEBUG_MODE
  printf("updateAlpha Function called with value %f\n", a);
#endif
  alpha = a;
  updateGUI();
}

void Color2GrayIt::updateMu(int m ){
#if DEBUG_MODE
  printf("updateMu Function called with value %d\n", m);
#endif
  r = m;
  updateGUI();
}
void Color2GrayIt::updateNumQColors(int n ){
#if DEBUG_MODE
  printf("NumQColors Function called with value %d\n", n);
#endif
  q_colors = n;
  updateGUI();
}




//CtrlWindow::CtrlWindow(){ CtrlW = make_window();}

//  gcc -c -g -D __APPLE__  -I/System/Library/Frameworks/Glut.framework/Versions/A/Headers/ -I/System/Library/Frameworks/OpenGL.framework/Versions/A/Headers/ -I/Users/gooch/Development/fltk-1.1.6 gui.cxx


