/*
 * This program is a modification of Peter Heckert's Unsharp masking plugin.
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 */


/*
 * This is a colour to greyscale converter for gimp-2.2.
 * The algorithm works first by mapping to Nayatani lightness in LUV
 * and then enhancing according to the chroma channel of the orignal image.
 *
 *
 * INSTALL with: $ gimptool-2.2  --install colour2grey-local-#.##.c 
 * (#.## is the version number)
 * (gimp-devel must be installed)
 *
 *
 * It will install in your user-plugin-folder 
 * To uninstall, simply delete it from  ~/.gimp-2.2/plugins
 * Or type gimptool --help and read and learn ;-)
 *
 * It will show up as "Colour2Grey Local" under the Filters -> Colors menu
 *
 *
 */ 

#define STANDALONE

#ifndef STANDALONE
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <gtk/gtk.h>
#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>

#ifdef STANDALONE

#define INIT_I18N voidproc
void voidproc(void){};
#define _(x) (x)
#define N_(x)  (x)

#else

#include "libgimp/stdplugins-intl.h"

#endif


#define PLUG_IN_VERSION "0.12"

#define SCALE_WIDTH   150
#define ENTRY_WIDTH     4

// This is for testing only!
#define float double
#define gfloat gdouble


/* Uncomment this line to get a rough estimate of how long the plug-in
 * takes to run.
 */

/*  #define TIMER  */


typedef struct
{
  gdouble  radius;
  gdouble  amount;
  gdouble  gamma;
  gboolean update_preview;
} C2gMaskParams;

typedef struct
{
  gboolean  run;
} Colour2GreyInterface;

/* local function prototypes */
static void      query (void);
static void      run   (const gchar      *name,
                        gint              nparams,
                        const GimpParam  *param,
                        gint             *nreturn_vals,
                        GimpParam       **return_vals);

static void      blur_line           (const gdouble  *ctable,
                                      const gdouble  *cmatrix,
                                      gint            cmatrix_length,
                                      const guchar   *src,
                                      guchar         *dest,
                                      gint            len,
                                      glong           bytes);

// Colour space conversions.
static void extract_lab (guchar  *src,
			 gint     bpp,
			 gint     numpix,
			 guchar *dst);

static void compose_lab (guchar *src,
			 gint     numpix,
			 gint    bytes,
			 guchar  *dst);

static void nayatani(guchar *src,
	     gint     bpp,
	     gint     numpix,
		     guchar *dst);

// Functions to run code.
static void
chromaunsharp (GimpPixelRgn *cPR,
               GimpPixelRgn *bcPR,
	       guchar  *g,
	       guchar  *bg,
	       gint     bpp,
	       gint x1,
	       gint y1,
	       gint     width,
	       gint     height,
	       gdouble     amount,
	       gdouble     gamma,
	       guchar   *dest);

static void      c2g_mask        (GimpDrawable   *drawable,
                                      gdouble         radius,
                                      gdouble         amount,
                                      gdouble         gamma);

static void      c2g_region      (GimpPixelRgn   *srcPTR,
                                      GimpPixelRgn   *dstPTR,
                                      gint            bytes,
                                      gdouble         radius,
                                      gdouble         amount,
                                      gdouble         gamma,
                                      gint            x1,
                                      gint            x2,
                                      gint            y1,
                                      gint            y2,
                                      gboolean        show_progress);

static gboolean  c2g_mask_dialog (GimpDrawable   *drawable);
static void      preview_update      (GimpPreview    *preview);


/* create a few globals, set default values */
static C2gMaskParams c2g_params =
  {
    5.0, /* default radius = 5 */
    0.0, /* default amount = .5 */
    1.0, /* gamma */
    TRUE /* default is to update the preview */
  };

/* Colour Globals */
static const double OFFSET_a    = -90.;
static const double OFFSET_b    =-130.;
static const double RANGE_a     = 190.;
static const double RANGE_b     = 230.;

// Reference white: D65
static const double ref_X = 0.950468;
static const double ref_Z = 1.08883;

/* Setting PLUG_IN_INFO */
GimpPlugInInfo PLUG_IN_INFO =
  {
    NULL,  /* init_proc  */
    NULL,  /* quit_proc  */
    query, /* query_proc */
    run,   /* run_proc   */
  };


// Gamma mapping conversion.
gfloat lut_gamma, inv_gamma;  
gfloat lut[256];  
  
static void 
lut_init(gfloat g)
{
 int i;
 
 for (i = 1; i< 255; i++) lut[i] = pow((gfloat) i/255.0, g);
 lut[0]= 0.0;
 lut[255] = 1.0;

 lut_gamma = g;
 inv_gamma = 1.0/g;
}


MAIN ()

static void
query (void)
{
  static GimpParamDef args[] =
    {
      { GIMP_PDB_INT32,    "run_mode",  "Interactive, non-interactive" },
      { GIMP_PDB_IMAGE,    "image",     "(unused)" },
      { GIMP_PDB_DRAWABLE, "drawable",  "Drawable to draw on" },
      { GIMP_PDB_FLOAT,    "radius",    "Radius of gaussian blur (in pixels > 1.0)" },
      { GIMP_PDB_FLOAT,    "amount",    "Strength of effect" },     
      { GIMP_PDB_FLOAT,    "gamma", "Gamma (param P in paper)" } ,            
    
    };

  gimp_install_procedure ("plug_in_colour2grey-local-"PLUG_IN_VERSION,
                          "A colour to greyscale converter",
                          "",
                          "",
			  "",
                          "",
			  N_("_Colour2Grey Local"),
                          "GRAY*, RGB*",
                          GIMP_PLUGIN,
                          G_N_ELEMENTS (args), 0,
                          args, NULL);

  gimp_plugin_menu_register ("plug_in_colour2grey-local-"PLUG_IN_VERSION,
                             "<Image>/Filters/Colors");
}

static void
run (const gchar      *name,
     gint              nparams,
     const GimpParam  *param,
     gint             *nreturn_vals,
     GimpParam       **return_vals)
{
  static GimpParam   values[1];
  GimpPDBStatusType  status = GIMP_PDB_SUCCESS;
  GimpDrawable      *drawable;
  GimpRunMode        run_mode;
#ifdef TIMER
  GTimer            *timer = g_timer_new ();
#endif

  run_mode = param[0].data.d_int32;

  *return_vals  = values;
  *nreturn_vals = 1;

  values[0].type          = GIMP_PDB_STATUS;
  values[0].data.d_status = status;

  INIT_I18N ();

  /*
   * Get drawable information...
   */
  drawable = gimp_drawable_get (param[2].data.d_drawable);
  gimp_tile_cache_ntiles (2 * (drawable->width / gimp_tile_width () + 1));

  switch (run_mode)
    {
    case GIMP_RUN_INTERACTIVE:
      gimp_get_data ("plug_in_colour2grey-local-"PLUG_IN_VERSION, &c2g_params);

      /* Reset default values show preview unmodified */
      /* initialize pixel regions and buffer */				// So basically, it seems that it shows the unmodified image.
      if (! c2g_mask_dialog (drawable))
        return;

      break;

    case GIMP_RUN_NONINTERACTIVE:
      if (nparams != 6)
        {
          status = GIMP_PDB_CALLING_ERROR;
        }
      else
        {
          c2g_params.radius = param[3].data.d_float;
          c2g_params.amount = param[4].data.d_float;
          c2g_params.gamma = param[5].data.d_float;          
          /* make sure there are legal values */
          if ((c2g_params.radius < 0.0) ||
              (c2g_params.amount < 0.0))
            status = GIMP_PDB_CALLING_ERROR;
        }
      break;

    case GIMP_RUN_WITH_LAST_VALS:
      gimp_get_data ("plug_in_colour2grey-local-"PLUG_IN_VERSION, &c2g_params);
      break;

    default:
      break;
    }

  if (status ==GIMP_PDB_SUCCESS)
    {
      drawable = gimp_drawable_get (param[2].data.d_drawable);

      /* here we go */
      c2g_mask (drawable, c2g_params.radius, c2g_params.amount, c2g_params.gamma);

      gimp_displays_flush ();

      /* set data for next use of filter */
      gimp_set_data ("plug_in_c2g_mask2-"PLUG_IN_VERSION, &c2g_params,
                     sizeof (C2gMaskParams));

      gimp_drawable_detach(drawable);
      values[0].data.d_status = status;
    }

#ifdef TIMER
  g_printerr ("%f seconds\n", g_timer_elapsed (timer, NULL));
  g_timer_destroy (timer);
#endif
}



static struct iir_param
{
 gdouble B,b1,b2,b3,b0,r,q;
 gdouble *p;
} iir; 

/*
static struct iir_param
{
gdouble B,b1,b2,b3,b0,r,q;
gdouble *p;
} iir; 
*/
static void iir_init(double r)
{
  
  iir.r = r;
  gdouble q;
  
  if ( r >= 2.5) q = 0.98711 * r - 0.96330;
  else q = 3.97156 - 4.14554 * sqrt(1.0-0.26891 * r);
  
  iir.q = q;
  iir.b0 = 1.57825 + ((0.422205 * q  + 1.4281) * q + 2.44413) *  q;
  iir.b1 = ((1.26661 * q +2.85619) * q + 2.44413) * q / iir.b0;
  iir.b2 = - ((1.26661*q +1.4281) * q * q ) / iir.b0;
  iir.b3 = 0.422205 * q * q * q / iir.b0;
  iir.B = 1.0 - (iir.b1 + iir.b2 + iir.b3);

}


static void iir_filter(gdouble *data, gint width)
/* 
 * Very fast gaussian blur with infinite impulse response filter
 * The row is blurred in forward direction and then in backward direction
 * So we achieve zero phase errors and symmetric impulse response
 * and good isotropy
 *
 * Theory for this filter can be found at:
 * <http://www.ph.tn.tudelft.nl/~lucas/publications/1995/SP95TYLV/SP95TYLV.pdf>
 * It is usable for radius downto 0.5. Lower radius must be done with the old
 * method. The old method also is very fast at low radius, so this doesnt matter 
 *
 * Double floating point precision is necessary for radius > 50, as my experiments
 * have shown. On my system (Duron, 1,2 GHz) the speed difference between double
 * and float is neglectable. 
 */

{

  gint w = c2g_params.radius;
  
  gdouble *const lp = data, *const rp = data + width-1+w;
  
  gint i;
    
  for (i=1; i<=w; i++) data[-i] = data[i]; /* mirror edges */
  for (i=1; i<=w; i++) data[i+width-1] = data[-i+width-1];
  
  {
  /* Hoping compiler will use optimal alternative, if not enough registers */
    register double d1,d2,d3;
    data = lp-w;
    d1=d2=d3=*data; 
    while (data <=  rp){
      *data *=  iir.B;
      *data +=  iir.b3 * d3;      
      *data +=  iir.b2 * (d3 = d2);    
      *data +=  iir.b1 * (d2 = d1); 
      d1 = *data++;
    } 
  
  
  }
  
  data--;   
  {
  /* Hoping compiler will use optimal alternative, if not enough registers */
    register double d1,d2,d3;
    d1=d2=d3=*data;
  
    while (data >=  lp){
      *data *=  iir.B;
      *data +=  iir.b3 * d3;      
      *data +=  iir.b2 * (d3 = d2);    
      *data +=  iir.b1 * (d2 = d1); 
      d1 = *data--;
    } 
  
  }
}

// -------------------- Gaussian Blur filter ---------------------
static void blur_line (const gdouble *ctable,
           const gdouble *cmatrix,
           gint           cmatrix_length,
           const guchar  *src,
           guchar        *dest,
           gint           len,   /* length of src and dest */
           glong          bytes) /* Bits per plane */	// ACtually Here is supposed to be how many bytes per pixel.
{


  gint    b = 0;
  gint    row;
  gint w = c2g_params.radius+10;

  // For each channel, dealing with double, so
  for (b = 0; b<bytes; b++){	// As many as the number of bytes per pixel
     gint idx;
	 //	idx = 0  ==> len
     for (row = b, idx=0 ; idx < len; row +=bytes, idx++) 
       iir.p[idx+w] = src[row];		// buffer

	 //	Later need to review
     iir_filter(iir.p+w, len);
     
     for (row =b, idx=0; idx < len; row +=bytes, idx++){
        gdouble value = iir.p[idx+w];     
        dest[row] = CLAMP( (gint) value, 0, 255); 
     }
  }
}

/***
 * 
 *  Functions for seperating the value information from the color 
 *  information. These are adapted from code in the compose and decompose
 *  plugins of Peter Kirchgessner.
 *
 */

const double Xn	= 0.951;
const double Yn	= 1.0;
const double Zn	= 1.089;

static void
chromaunsharp (GimpPixelRgn *cPR,
               GimpPixelRgn *bcPR,
	       guchar  *g,
	       guchar  *bg,
	       gint     bpp,
	       gint x1,
	       gint y1,
	       gint width,
	       gint height,
	       gdouble amount,
	       gdouble gamma,
	       guchar   *dest)
{

  register guchar *lab_dest = dest;

  guchar *c = g_new(guchar, width*height*bpp);
  guchar *bc = g_new(guchar, width*height*bpp);
  gimp_pixel_rgn_get_rect (cPR, c, x1, y1, width, height);	// i guess it gets the original image..
  gimp_pixel_rgn_get_rect (bcPR, bc, x1, y1, width, height); // i guess it gets the blurred colour image.
  
  guchar *ptr[4];
  ptr[0] = c;	// original image data pointer.
  ptr[1] = bc;	// original blurred image data pointer.
  ptr[2] = g;	// nayatani Grey in Lab format(?) - lab? not sure...
  ptr[3] = bg;	// nayatani blurred grey in Lab format- lab? not sure..

  register gint count = width*height, offset = bpp; // bpp = bytes per pixel

  gint i;
  gdouble red, green, blue;
  gdouble x, y, z;
  gdouble tl, l[4], a[4], b[4];
  gdouble tx, ty, tz;
  gdouble ftx, fty, ftz;
  gdouble lambda, de, hg, deg;
  gdouble sixteenth = 16.0 / 116.0;

  while (count-- > 0)
    {
    for ( i=0; i < 4; i++ ) {
      // C to RGB
      red   = *ptr[i] / 255.0;	/// ptr[0] == C ptr[1] == bc ptr[2] == g ptr[3] == bg;
      green = *ptr[i] / 255.0;
      blue  = *ptr[i] / 255.0;

      //apply gamma function to convert to sRGB
      if ( red <=0.03928 ) red = red/12.92;
      else red = pow((red+0.055)/1.055,2.4);
      if ( green <=0.03928 ) green = green/12.92;
      else green = pow((green+0.055)/1.055,2.4);
      if ( blue <=0.03928 ) blue = blue/12.92;
      else blue = pow((blue+0.055)/1.055,2.4);

      x = 0.431 * red + 0.342 * green + 0.178 * blue;
      y = 0.222 * red + 0.707 * green + 0.071 * blue;
      z = 0.020 * red + 0.130 * green + 0.939 * blue;

      if (( ty = y/Yn ) > 0.008856)
        {
          tl   = 116.0 * cbrt( ty ) - 16.0;
          fty = cbrt( ty );
        }
      else
        {
          tl   = 903.3 * ty;
          fty = 7.78*ty + sixteenth;
        }

      ftx = ((tx = x/Xn) > 0.008856) ? cbrt (tx) : 7.78 * tx + sixteenth;
      ftz = ((tz = z/Zn) > 0.008856) ? cbrt (tz) : 7.78 * tz + sixteenth;

      l[i] = tl;
      a[i] = ftx - fty;
      b[i] = fty - ftz;

      ptr[i] += offset;
    }  // end for

	// now lab.

    // Calculate lambda
    hg = l[2]-l[3];			// G bandpass == h(G)
    de = pow(pow(l[0]-l[1],2)+pow(a[0]-a[1],2)+pow(b[0]-b[1],2),0.5);	// de = h(I)
    deg = pow(pow(hg,2),0.5);
    if ( deg == 0 ) 
      lambda = de/0.00001;
     else
       lambda = de/deg;

    lambda = pow(lambda,gamma);
    tl = l[2] + (amount*lambda*hg);			// I[2] == nayatani G.

    // Keep the old chromatic greyimage's channels.
    lab_dest[0] = CLAMP( (gint)( tl*2.5599) , 0, 255);
    lab_dest[1] = (guchar) (128.0 + a[2] * 635 );	// I think a and b have been set to '128' in function nayatani()
    lab_dest[2] = (guchar) (128.0 + b[2] * 254 ); 
    lab_dest += offset;
    
  }

  g_free (c);
  g_free (bc);
}


static void
extract_lab (guchar  *src,
	     gint     bpp,
	     gint     numpix,
	     guchar *dst)
{

  register guchar *rgb_src = src;
  register guchar *lab_dst = dst;
  register gint count = numpix, offset = bpp;

  gdouble red, green, blue;
  gdouble x, y, z;
  gdouble l; /*, a, b; */
  gdouble tx, ty, tz;
  gdouble ftx, fty, ftz;

  gdouble sixteenth = 16.0 / 116.0;

  while (count-- > 0)
    {
      red   = rgb_src[0] / 255.0;
      green = rgb_src[1] / 255.0;
      blue  = rgb_src[2] / 255.0;

      //apply gamma function to convert to sRGB
      if ( red <=0.03928 ) red = red/12.92;
      else red = pow((red+0.055)/1.055,2.4);
      if ( green <=0.03928 ) green = green/12.92;
      else green = pow((green+0.055)/1.055,2.4);
      if ( blue <=0.03928 ) blue = blue/12.92;
      else blue = pow((blue+0.055)/1.055,2.4);

      x = 0.431 * red + 0.342 * green + 0.178 * blue;
      y = 0.222 * red + 0.707 * green + 0.071 * blue;
      z = 0.020 * red + 0.130 * green + 0.939 * blue;

      if (( ty = y/Yn ) > 0.008856)
        {
          l   = 116.0 * cbrt( ty ) - 16.0;
          fty = cbrt( ty );
        }
      else
        {
          l   = 903.3 * ty;
          fty = 7.78*ty + sixteenth;
        }

      ftx = ((tx = x/Xn) > 0.008856) ? cbrt (tx) : 7.78 * tx + sixteenth;
      ftz = ((tz = z/Zn) > 0.008856) ? cbrt (tz) : 7.78 * tz + sixteenth;

      lab_dst[0] = (guchar) (l * 2.5599);
      lab_dst[1] = (guchar) (128.0 + (ftx - fty) * 635 );
      lab_dst[2] = (guchar) (128.0 + (fty - ftz) * 254 );

      rgb_src += offset;
      lab_dst += offset;
    }
}

/*
basically convert to sRGB
*/
static void
compose_lab (guchar *src,
             gint    numpix,
	     gint    bytes,
             guchar  *dst)
{
  register guchar *lab_src = src;
  register guchar *rgb_dst = dst;

  register gint count = numpix;

  gdouble red, green, blue;
  gdouble x, y, z;
  gdouble l, a, b;

  gdouble p, yyn;
  gdouble ha, hb, sqyyn;
  gdouble g=1/2.4;

  while (count-- > 0)
    {
      l = lab_src[0] / 2.550;
      a = ( lab_src[1] - 128.0 ) / 1.27;
      b = ( lab_src[2] - 128.0 ) / 1.27;

      p = (l + 16.) / 116.;
      yyn = p*p*p;

      if (yyn > 0.008856)
        {
          y = Yn * yyn;
          ha = (p + a/500.);
          x = Xn * ha*ha*ha;
          hb = (p - b/200.);
          z = Zn * hb*hb*hb;
        }
      else
        {
          y = Yn * l/903.3;
          sqyyn = pow(l/903.3,1./3.);
          ha = a/500./7.787 + sqyyn;
          x = Xn * ha*ha*ha;
          hb = sqyyn - b/200./7.787;
          z = Zn * hb*hb*hb;
        };

      red   =  3.063 * x - 1.393 * y - 0.476 * z;
      green = -0.969 * x + 1.876 * y + 0.042 * z;
      blue  =  0.068 * x - 0.229 * y + 1.069 * z;

      if ( red <=0.00304 ) red = red*12.92;
      else red = pow(red*1.055,g)-0.055;
      if ( green <=0.00304 ) green = green*12.92;
      else green = pow(green*1.055,g)-0.055;
      if ( blue <=0.00304 ) blue = blue*12.92;
      else blue = pow(blue*1.055,g)-0.055;

      red   = ( red   > 0 ) ? red   : 0;
      green = ( green > 0 ) ? green : 0;
      blue  = ( blue  > 0 ) ? blue  : 0;

      red   = ( red   < 1.0 ) ? red   : 1.0;
      green = ( green < 1.0 ) ? green : 1.0;
      blue  = ( blue  < 1.0 ) ? blue  : 1.0;

      
      rgb_dst[0] = (guchar) ( red   * 255.999 );
      rgb_dst[1] = (guchar) ( green * 255.999 );
      rgb_dst[2] = (guchar) ( blue  * 255.999 );

      rgb_dst += bytes;
      lab_src += bytes;
    }
}

static void nayatani(guchar *src,
	     gint     bpp,
	     gint     numpix,
	     guchar *dst)
{

  register guchar *lab_src = src;
  register guchar *lab_dst = dst;
  register gint count = numpix;

  gdouble x, y, z;
  gdouble l, a, b;

  gdouble p, yyn;
  gdouble ha, hb, sqyyn;

  gdouble u, v, kbr, hue, qhue, suv, gamma;
  gdouble adaptlum = 20.0;
  gdouble whiteu = 0.20917;
  gdouble whitev = 0.48810;
  // Rwhite xyz:  96.420  100.000   82.490

  while (count-- > 0)
    {
      l = lab_src[0] / 2.550;
      a = ( lab_src[1] - 128.0 ) / 1.27;
      b = ( lab_src[2] - 128.0 ) / 1.27;

      p = (l + 16.) / 116.;
      yyn = p*p*p;

      if (yyn > 0.008856)
        {
          y = Yn * yyn;
          ha = (p + a/500.);
          x = Xn * ha*ha*ha;
          hb = (p - b/200.);
          z = Zn * hb*hb*hb;
        }
      else
        {
          y = Yn * l/903.3;
          sqyyn = pow(l/903.3,1./3.);
          ha = a/500./7.787 + sqyyn;
          x = Xn * ha*ha*ha;
          hb = sqyyn - b/200./7.787;
          z = Zn * hb*hb*hb;
        };

      u=4*x/(x+15*y+3*z)-whiteu;
      v=9*y/(x+15*y+3*z)-whitev;
      hue=atan2(v,u);
      qhue = -0.01585-0.03016*cos(hue)-0.04556*cos(2*hue)- 0.02667*cos(3*hue)-0.00295*cos(4*hue)+0.14592*sin(hue)+0.05084*sin(2*hue)-0.01900*sin(3*hue)-0.00764*sin(4*hue);

      // Adapting luminance dependency
      kbr = 0.2717*pow(6.469+6.362*adaptlum,0.4495)/pow(6.469+adaptlum,0.4495);

      // Saturation
      suv=13*pow(pow(u,2) + pow(v,2),0.5);

      // VAC from equation 7
      gamma = 1 + (-0.1340*qhue + 0.0872*kbr)*suv;
      //gamma = 1 + (-0.8660*qhue + 0.0872*kbr)*suv; // VCC

      lab_dst[0] = CLAMP( (gint)(gamma*l * 2.5599) , 0, 255);
      lab_dst[1] = (guchar) 128.0;
      lab_dst[2] = (guchar) 128.0;

      lab_src += bpp;
      lab_dst += bpp; // bpp == the number of bytes per pixel... so it is just moving on to the next pixel.
    }
}

//-------------------------------------------------------------------------
static void
c2g_mask (GimpDrawable *drawable,
              gdouble       radius,
              gdouble       amount,
              gdouble       gamma)
{

  GimpPixelRgn srcPR, destPR;
  gint         x1, y1, x2, y2;

  /* initialize pixel regions */
  gimp_pixel_rgn_init (&srcPR, drawable,
                       0, 0, drawable->width, drawable->height, FALSE, FALSE);
  gimp_pixel_rgn_init (&destPR, drawable,
                       0, 0, drawable->width, drawable->height, TRUE, TRUE);

  /* Get the input */ //  Find the bounding box of the current selection in relation to the specified drawable. 
  gimp_drawable_mask_bounds (drawable->drawable_id, &x1, &y1, &x2, &y2);

  c2g_region (&srcPR, &destPR, drawable->bpp,
                  radius, amount, gamma,
                  x1, x2, y1, y2,
                  TRUE);

  gimp_drawable_flush (drawable);
  gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
  gimp_drawable_update (drawable->drawable_id, x1, y1, x2 - x1, y2 - y1);
}

//------------------------------------------------------------------------
/* Perform an c2g mask on the region, given a source region, dest.
 * region, width and height of the regions, and corner coordinates of
 * a subregion to act upon.  Everything outside the subregion is unaffected.
 */

static void
c2g_region (GimpPixelRgn *srcPR,
                GimpPixelRgn *destPR,
                gint          bytes, /* Bytes per pixel */
                gdouble       radius,
                gdouble       amount_p,
                gdouble       gamma_p,
                gint          x1,    /* Corners of subregion */
                gint          x2,
                gint          y1,
                gint          y2,
                gboolean      show_progress)
{
  guchar  *src;
  guchar  *dest;
  guchar *bigsrc;
  guchar *bigdest;
  guchar *tmpdest;
 
  gint     width   = x2 - x1;
  gint     height  = y2 - y1;
  gdouble *cmatrix = NULL;
  gint     cmatrix_length;
  gdouble *ctable;
  gint     row, col, idx;
  gint w = c2g_params.radius+10;

  gdouble amount = c2g_params.amount;
  gdouble gamma = c2g_params.gamma;
    
  if (show_progress)
    gimp_progress_init (_("Blurring..."));
  /* gimp_progress_init
  Initializes the progress bar for the current plug-in. 

  Initializes the progress bar for the current plug-in. It is only valid to call this procedure from a plug-in.

  message :
  Message to use in the progress dialog.  
  Returns :
  TRUE on success.  
  */
  
  iir_init(c2g_params.radius);				// I don't really know what it does.

  /* allocate buffers */
  src  = g_new (guchar, MAX (width, height) * bytes /* bytes per pixel*/);
  dest = g_new (guchar, MAX (width, height) * bytes);
  iir.p = g_new(gdouble, MAX (width, height)+2*w);

 
  bigsrc = g_new(guchar, width * height * bytes);
  bigdest = g_new(guchar, width * height * bytes);
  tmpdest = g_new(guchar, width * height * bytes);

  if (show_progress)
    gimp_progress_init (_("Colour converting..."));
  
  // 1. Calculate Nayatani Grey and Blur in LAB
  gimp_pixel_rgn_get_rect (srcPR, bigsrc, x1, y1, width, height);
  extract_lab(bigsrc, bytes, width * height, bigdest);
  nayatani(bigdest,bytes,width * height, bigdest);	// From this point, bigdest is I guess a Grey image L + a(= 128) + b(= 128)
  gimp_pixel_rgn_set_rect (destPR, bigdest, x1, y1, width, height);	// bigdest == destPR == nayatani grey

  // 2. Make a blur of Grey LAB
  //	Basically it blurs each row every loop.
  for (row = 0, idx=0; row < height; row++, idx+=width)
    {
      gimp_pixel_rgn_get_row (destPR, src, x1, y1 + row, width);
      blur_line (ctable, cmatrix, cmatrix_length, src, dest, width, bytes);   
      gimp_pixel_rgn_set_row (destPR, dest, x1, y1 + row, width);
    }

	//	Basically it blurs each column every loop.
  for (col = 0; col < width; col++)
    {
      gimp_pixel_rgn_get_col (destPR, src, x1 + col, y1, height);
      blur_line (ctable, cmatrix, cmatrix_length, src, dest, height, bytes);
      gimp_pixel_rgn_set_col (destPR, dest, x1 + col, y1, height);

	  //	no big deal, it just updates the progress bar which has been created from c2g_mask_dialog function()
      if (show_progress && col % 8 == 0)
        gimp_progress_update ((gdouble) col / (3 * width) + 0.33);
    }

  // 3. Convert grey and blur back to RGB. at this point, destPR is a blurred grey image
	//	destPR is the GimpPixelRgn type.
  compose_lab(bigdest, width * height, bytes, bigdest); // bigdest=greyRGB
  gimp_pixel_rgn_get_rect (destPR, bigsrc, x1, y1, width, height);
  compose_lab(bigsrc, width * height, bytes, bigsrc); // bigsrc= out of bigdest GimpPixelRgn

  // 4. Blur Colour RGB and write into destPR
  gimp_pixel_rgn_get_rect (srcPR, tmpdest, x1, y1, width, height);
  gimp_pixel_rgn_set_rect (destPR, tmpdest, x1, y1, width, height);

  // I guess it blurs the grey image in RGB format -> NO.,....
  for (row = 0, idx=0; row < height; row++, idx+=width)
    {
      gimp_pixel_rgn_get_row (destPR, src, x1, y1 + row, width);
      blur_line (ctable, cmatrix, cmatrix_length, src, dest, width, bytes);   
      gimp_pixel_rgn_set_row (destPR, dest, x1, y1 + row, width);
    }

  for (col = 0; col < width; col++)
    {
      gimp_pixel_rgn_get_col (destPR, src, x1 + col, y1, height);
      blur_line (ctable, cmatrix, cmatrix_length, src, dest, height, bytes);
      gimp_pixel_rgn_set_col (destPR, dest, x1 + col, y1, height);

      if (show_progress && col % 8 == 0)
        gimp_progress_update ((gdouble) col / (3 * width) + 0.33);
    }

  // destPR = blur colour RGB
  chromaunsharp( srcPR, destPR, bigdest, bigsrc, bytes, x1, y1, width, height, amount, gamma, tmpdest );
  // tmpdest is the result of chromaunsharp
  compose_lab(tmpdest, width * height, bytes, tmpdest); // tmpdest has unsharp
  //	tmpdest is G' now.
  gimp_pixel_rgn_set_rect (destPR, tmpdest, x1, y1, width, height);
  

  if (show_progress)
    gimp_progress_update (0.0);

  g_free (bigsrc);
  g_free (bigdest);
  g_free (tmpdest);
  g_free (iir.p);
  g_free (dest);
  g_free (src);
  
}

/* generates a 1-D convolution matrix to be used for each pass of
 * a two-pass gaussian blur.  Returns the length of the matrix.
 */

static gboolean
c2g_mask_dialog (GimpDrawable *drawable)
{
  GtkWidget *dialog;
  GtkWidget *main_vbox;
  GtkWidget *preview;
  GtkWidget *table;
  GtkObject *adj;
  gboolean   run;

  gimp_ui_init ("c2g-"PLUG_IN_VERSION, TRUE);

  dialog = gimp_dialog_new (_("C2G Global/Local "PLUG_IN_VERSION), "c2g2-"PLUG_IN_VERSION,
                            NULL, 0,
                            gimp_standard_help_func, "plug-in-c2g-mask2-"PLUG_IN_VERSION,

                            GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
                            GTK_STOCK_OK,     GTK_RESPONSE_OK,

                            NULL);

  main_vbox = gtk_vbox_new (FALSE, 4);
  gtk_container_set_border_width (GTK_CONTAINER (main_vbox), 4);
  gtk_container_add (GTK_CONTAINER (GTK_DIALOG (dialog)->vbox), main_vbox);
  gtk_widget_show (main_vbox);

  preview = gimp_drawable_preview_new (drawable,
                                       &c2g_params.update_preview);
  gtk_box_pack_start (GTK_BOX (main_vbox), preview, TRUE, TRUE, 0);
  gtk_widget_show (preview);

  g_signal_connect (preview, "invalidated",
                    G_CALLBACK (preview_update),
                    NULL);

  table = gtk_table_new (3, 3, FALSE);
  gtk_table_set_col_spacings (GTK_TABLE (table), 6);
  gtk_table_set_row_spacings (GTK_TABLE (table), 6);
  gtk_box_pack_start (GTK_BOX (main_vbox), table, FALSE, FALSE, 0);
  gtk_widget_show (table);

  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, 0,
                              _("_Radius:"), SCALE_WIDTH, ENTRY_WIDTH,
                              c2g_params.radius, 0.25, 500.0, 0.1, 1.0, 1,
                              TRUE, 0, 0,
                              NULL, NULL);
  gimp_scale_entry_set_logarithmic(adj,TRUE);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &c2g_params.radius);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);

  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, 1,
                              _("_Amount +:"), SCALE_WIDTH, ENTRY_WIDTH,
                              c2g_params.amount, 0.0, 5.0, 0.01, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &c2g_params.amount);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);

  adj = gimp_scale_entry_new (GTK_TABLE (table), 0, 2,
                              "_Gamma:", SCALE_WIDTH, ENTRY_WIDTH,
                              c2g_params.gamma,
                              0.1, 5.0, 0.01, 0.1, 2,
                              TRUE, 0, 0,
                              NULL, NULL);

  g_signal_connect (adj, "value_changed",
                    G_CALLBACK (gimp_double_adjustment_update),
                    &c2g_params.gamma);
  g_signal_connect_swapped (adj, "value_changed",
                            G_CALLBACK (gimp_preview_invalidate),
                            preview);

//-----------------------------------------------------------------------------  
  gtk_widget_show (dialog);

  run = (gimp_dialog_run (GIMP_DIALOG (dialog)) == GTK_RESPONSE_OK);

  gtk_widget_destroy (dialog);

  return run;
}

static void
preview_update (GimpPreview *preview)
{
  GimpDrawable *drawable;
  gint          x1, x2;
  gint          y1, y2;
  gint          x, y;
  gint          width, height;
  gint          border;
  GimpPixelRgn  srcPR;
  GimpPixelRgn  destPR;

  drawable =
    gimp_drawable_preview_get_drawable (GIMP_DRAWABLE_PREVIEW (preview));

  gimp_pixel_rgn_init (&srcPR, drawable,
                       0, 0, drawable->width, drawable->height, FALSE, FALSE);
  gimp_pixel_rgn_init (&destPR, drawable,
                       0, 0, drawable->width, drawable->height, TRUE, TRUE);

  gimp_preview_get_position (preview, &x, &y);
  gimp_preview_get_size (preview, &width, &height);

  /* enlarge the region to avoid artefacts at the edges of the preview */
  border = 2.0 * c2g_params.radius + 0.5;
  
  if (border > width/2) border = width/2; //Speed up preview 
  x1 = MAX (0, x - border);
  y1 = MAX (0, y - border);
  x2 = MIN (x + width  + border, drawable->width);
  y2 = MIN (y + height + border, drawable->height);

  c2g_region (&srcPR, &destPR, drawable->bpp,
                  c2g_params.radius, c2g_params.amount, c2g_params.gamma,
                  x1, x2, y1, y2,
                  FALSE);

  gimp_pixel_rgn_init (&destPR, drawable, x, y, width, height, FALSE, TRUE);
  gimp_drawable_preview_draw_region (GIMP_DRAWABLE_PREVIEW (preview), &destPR);
}
