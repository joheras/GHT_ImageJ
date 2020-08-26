package unirioja.GHT;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import ij.ImagePlus;
import net.imagej.Dataset;
import net.imagej.ImageJ;
import net.imagej.display.ImageDisplay;
import net.imagej.display.ImageDisplayService;
import net.imagej.ops.OpService;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.IntegerType;

@Plugin(type = Command.class, headless = true, menuPath = "Image>Adjust>Threshold GHT")
public class GHT implements Command {

	@Parameter(min = "0")
	private float nu = 0;

	@Parameter(min = "0")
	private float tau = 0;

	@Parameter(min = "0")
	private float kappa = 0;

	@Parameter(min = "0", max = "1")
	private float omega;

	/*
	 * This {@code @Parameter} is for Image Display service. The context will
	 * provide it automatically when this command is created.
	 */
	@Parameter
	private ImageDisplayService imageDisplayService;

	/*
	 * This {@code @Parameter} is for ImageJ Ops service. The context will provide
	 * it automatically when this command is created.
	 */
	@Parameter
	private OpService opService;

	@Parameter(type = ItemIO.INPUT)
	private ImageDisplay displayIn;

	private static double[] csum(long[] histo) {
		double[] res = new double[histo.length - 1];

		res[0] = histo[0];
		for (int i = 1; i < res.length; i++) {
			res[i] = histo[i] + res[i - 1];
		}

		return res;
	}
	
	private static double[] csum(double[] histo) {
		double[] res = new double[histo.length - 1];

		res[0] = histo[0];
		for (int i = 1; i < res.length; i++) {
			res[i] = histo[i] + res[i - 1];
		}

		return res;
	}

	private static double[] dsum(long[] histo) {
		double[] res = new double[histo.length];
		double[] res1 = new double[histo.length - 1];

		res[0] = histo[histo.length - 1];
		for (int i = 1; i < res.length; i++) {
			res[i] = histo[histo.length - 1 - i] + res[i - 1];
		}

		for (int i = 0; i < res1.length; i++) {
			res1[i] = res[res1.length - i - 1];
		}

		return res1;
	}
	
	private static double[] dsum(double[] histo) {
		double[] res = new double[histo.length];
		double[] res1 = new double[histo.length - 1];

		res[0] = histo[histo.length - 1];
		for (int i = 1; i < res.length; i++) {
			res[i] = histo[histo.length - 1 - i] + res[i - 1];
		}

		for (int i = 0; i < res1.length; i++) {
			res1[i] = res[res1.length - i - 1];
		}

		return res1;
	}

	private static double[] clip(double[] histo) {
		double[] res = new double[histo.length];
		for (int i = 0; i < res.length; i++) {
			res[i] = Math.max(histo[i],(long)Math.pow(10,-30));
		}
		return res;
	}

	private static double[] arange(long[] histo) {

		double[] res = new double[histo.length];
		for (int i = 0; i < res.length; i++) {
			res[i] = i;
		}
		return res;

	}
	
	private static int ght_threshold(long[] histogram,float nu,float tau,float kappa,float omega) {
		double[] w0 = clip(csum(histogram));
		double[] w1 = clip(dsum(histogram));

		double[] w0w1 = new double[w0.length];

		for (int i = 0; i < w0.length; i++) {
			
			w0w1[i] = w0[i]+ w1[i];
		}

		double[] p0 = new double[w0.length];
		double[] p1 = new double[w0.length];

		for (int i = 0; i < w0.length; i++) {
			
			p0[i] = w0[i]/ w0w1[i];
			p1[i] = w1[i]/w0w1[i];
		}

		double[] x = arange(histogram);

		double[] nx = new double[histogram.length];
		double[] nxx = new double[histogram.length];
		for (int i = 0; i < histogram.length; i++) {
			nx[i] = histogram[i]* x[i];
			nxx[i] = histogram[i]* x[i] * x[i];
		}

		double[] mu0 = csum(nx);
		double[] mu1 = dsum(nx);
		for (int i = 0; i < w0.length; i++) {
			mu0[i] = mu0[i]/ w0[i];
			mu1[i] = mu1[i]/w1[i];
		}

		double[] d0 = csum(nxx);
		double[] d1 = dsum(nxx);

		for (int i = 0; i < w0.length; i++) {
			d0[i] -= w0[i] * mu0[i] * mu0[i];
			d1[i] -= w1[i] * mu1[i] * mu1[i];
		}

		double[] v0 = new double[p0.length];
		double[] v1 = new double[p0.length];
		double[] f0 = new double[p0.length];
		double[] f1 = new double[p0.length];
		double[] f0f1 = new double[p0.length];

		for (int i = 0; i < p0.length; i++) {
			v0[i] = Math.max((double)(p0[i] *  nu *  tau * tau + d0[i]) / (p0[i] * nu + w0[i]),Math.pow(10,-30));
			v1[i] = Math.max((p1[i] * nu *  tau *  tau + d1[i]) / (p1[i] *  nu + w1[i]),Math.pow(10,-30));
			
			
			
			f0[i] = -d0[i] / v0[i] - w0[i] * Math.log(v0[i])
					+ 2 * (w0[i] +  kappa *  omega) * Math.log(w0[i]);
			f1[i] = -d1[i] / v1[i] - w1[i] *  Math.log(v1[i])
					+ 2 * (w1[i] + kappa *  (1 - omega)) *  Math.log(w1[i]);
			f0f1[i] = f0[i] + f1[i];
			if(Double.isNaN(f0f1[i])) {f0f1[i]=0;}
		}

		int sumMax = 0;
		int nmax = 1;
		double max = f0f1[0];

		for (int i = 1; i < f0f1.length; i++) {
			if (f0f1[i] == max) {
				nmax++;
				sumMax += i;
			}
			if (f0f1[i] > max) {
				max = f0f1[i];
				nmax = 1;
				sumMax = i;
			}
		}

		return (sumMax / nmax);
	}
	
	

	public void run() {
		final Dataset input = imageDisplayService.getActiveDataset(displayIn);
		Img<IntegerType> image = (Img<IntegerType>) input.getImgPlus();
		Histogram1d<IntegerType> hist = opService.image().histogram(image);
		long[] histogram = hist.toLongArray();
		System.out.println(ght_threshold(histogram,0,0,0,0));


		

	}



}
