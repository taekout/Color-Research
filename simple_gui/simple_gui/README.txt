README.txt for Color2Gray simple GUI
by Amy Gooch

The underlying code is the same as the simple source code.
I've built a FLTK interface, which allows you to control the variables
associated with our method, including the new speed-up variables
(which uses color quanitzation).

--------> DISCLAIMER:
This is development code, not comercial code.. it may be full of bugs,
(but they are cute bugs.. you know, the ones with spots). At the
least, it may be a little bit of frustrating.  

This code is Copyright (c) 2005 Amy Gooch and  Sven Olsen
Patent-pending.
Provided that we are not held responsible for any resulting damages,
you are hereby given permission to modify and redistribute it as you
like.

--------> REQUIREMENTS:
	fltk :  http://www.fltk.org/
	      (modify the makefile to include your fltk path)
	image magick commandline tools: 
	      http://www.imagemagick.org/script/index.php

Note: This was originally built on a Mac, using OS X.3, so the
makefile is overly complicated and not very clean.

--------> COMPILING:

On a mac, you need to type:
   make; make macosx

(the extra step associates the picky resource files, otherwise you
will notice on your mac that the app is running, but only in the
background).

On other machines:
   well... you know what to do....

--------> RUNNING:

./color2gray 
	     and then just use the interface
OR
you can give it the same arguments as the simple_source program:
    ./color2gray  Sunrise.ppm -theta 45 -alpha 10 -mu 0 -q 64
 

--------> THE INTERFACE:

First thing to do is to press
      LoadSource

This will pop up a browser window & show you all of the ppm files
available
By default the program will load up cb45sm.ppm included in this directory

Next you can press
     Run Color2Gray
After some time (depending upon the size of your images), it will pop
     up a window with the original color image, default grayscale
     image, new color2gray image, and new color2gray+chrominance
     image, using the parameters set in the yellow window.

Subsequent presses to 
	   Run Color2Gray
     will recalculate the images, and refresh the window based upon changes made to the values in
	   the yellow box.

Parameters:
----------
theta: controls how chrominance differences are mapped to positive or
	negative differences in luminance space (DEFAULT = 45)
alpha: controls the scale of the chrominance offset (DEFAULT = 10)
mu: radius of pixels considered as neighbors (DEFAULT = 0 = Full
	neighbor hood)
quantize (check box): Use speed-up method, requires mu=0, which quantizes
	the color image to use "Num Colors".
	See:
	http://www.cs.northwestern.edu/~sco590/c2g_notes/notes.pdf
	for details.

It is suggested that you test on images about 200x200, especially if
you are not using quanization method


If the theta parameter is not making a lot of sense, we provide a
visualization of the input image scaled down to 200x200 or smaller
image, over various values of the parameter theta.
Currently the code assigns the variables to the values:
	  mu = 0
	  alpha = 10
	  Quantize = true
	  num_colors = 30
And the user interface will ~not~ change these values.

All of these images are saved out in the program directory with the formula:
    sprintf(outname, "c2g_theta%.1f_a%.1f_mu%s_q%d.%d_%s",
                             theta*r2d,alpha,mu,quantize,q_colors,fname)
			    
   where fname is the original filename.

  Color2Gray+Chrominance images are also saved out, but start with
  "c2gC_" instead of "c2g_"


--------> List of To Dos:

Control over image size for "Run Color2Gray"

Smarter placement of result images

Control over values for theta wheel
Control over image size for theta wheel

