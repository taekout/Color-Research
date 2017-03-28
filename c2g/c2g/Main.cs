// c2g -- Color to gray conversion
// (C)Copyright 2005 by martin.faust@e56.de
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the
// Free Software Foundation, Inc., 59 Temple Place - Suite 330,
// Boston, MA  02111-1307, USA.
using System;
using System.Drawing;

namespace c2g {
	class MainClass {
		public static void Main(string[] args) {
			int x, y;
			float min, max, mean;
			Color pixel;
			float saturation, brightness, hue;
			float[,] gray;
			bool CorrectBrightness;
			
			Console.WriteLine("c2g by martin.faust@e56.de");
			if (args == null || args.Length != 3) {
				Console.WriteLine("usage: c2g <input rgb name> <output gray name> <CorrectBrightness: true|false>");
				return;
			}
			
			CorrectBrightness = args[2].Equals("true");

			Console.WriteLine("Reading {0}", args[0]);
			Bitmap bm = Image.FromFile(args[0]) as Bitmap;
			
			Console.WriteLine("Processing {0}x{1}", bm.Width, bm.Height);
			gray = new float[bm.Width,bm.Height];
			min = 100.0f;
			max = -100.0f;
			mean = 0.0f;
			for(y=0; y<bm.Height; y++) {
				for(x=0; x<bm.Width; x++) {
					pixel = bm.GetPixel(x, y);
					saturation = pixel.GetSaturation();
					brightness = pixel.GetBrightness();
					hue = pixel.GetHue();
					
					if (saturation == 0.0f)
						gray[x,y] = 1.5f * brightness;
					else
						gray[x, y] = brightness +  brightness*saturation;
					
					if (gray[x, y] < min) min = gray[x, y];
					if (gray[x, y] > max) max = gray[x, y];
					mean += gray[x, y];
				}
			}
			
			mean /= (float) (bm.Height * bm.Width);
			Console.WriteLine("Grey values: min {1} mean {2} max {3}", args[1], min, mean, max);
			min = 0.0f; max = (mean + max) * 0.5f;
	
			Console.WriteLine("Creating gray scale image");			
			for(y=0; y<bm.Height; y++) {
				for(x=0; x<bm.Width; x++) {
					if (CorrectBrightness)
						brightness = 0.9f * 255.0f * (gray[x, y] - min) / (max - min);
					else
						brightness = 255.0f * (gray[x, y] - min) / (max - min);
					if (brightness > 255.0f) brightness = 255.0f;
					if (brightness < 0.0f) brightness = 0.0f;
					bm.SetPixel(x, y, Color.FromArgb(0, (int) brightness, (int) brightness, (int) brightness));
				}
			}
			
			Console.WriteLine("Saving as {0}", args[1]);
			bm.Save(args[1]);
		}
	}
}
