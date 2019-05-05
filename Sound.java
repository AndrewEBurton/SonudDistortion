/**
* Title: Audio Analysis Using a FFT
* Author: Andrew Burton
* Current version written: May 2016
* Description: FT of a pure tone 
*/

// Import packages
import java.awt.Color;
import java.awt.Font;
import java.io.File;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.math.plot.*;
import org.math.plot.plotObjects.*;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

import org.math.plot.Plot2DPanel;
import org.math.plot.plotObjects.BaseLabel;

public class SoundDistort {

	// Start main method
	public static void main(String[] args) {

		// Input file name
		String wavFileName = "./Audio/C5Cello.wav";		

		// Initialize things outside the try/catch	
		double dt = 0.0;
		double fNyq = 0.0;
		double[] time = null;
		double[] signal = null;
        double[] signalDist = null;
        long sampleRate = 0;
        int numFrames = 0;
        
        
		//////////////////////////////////////		
		// 1. Reading WAV file
		//////////////////////////////////////

		try {		
			// Open the WAV file
			WavFile wavFile = WavFile.openWavFile(new File(wavFileName));

			// Display information about the WAV file	
			wavFile.display();
		
			// Calculate the time step, Nyquist frequency, Sample Rate, and Number of frames
			dt = (double)(1./wavFile.getSampleRate());
			fNyq = (double)(wavFile.getSampleRate() / 2.);
            sampleRate = wavFile.getSampleRate();
            numFrames = (int)wavFile.getNumFrames();

			// Set the size of the arrays that hold the signal and time values
            signal = new double[(int)wavFile.getNumFrames()];	
            signalDist = new double[(int)wavFile.getNumFrames()];
			time = new double[(int)wavFile.getNumFrames()];
			
			// Read the WAV file
			wavFile.readFrames(signal,signal.length);

			// Create and fill the time array
			for (int n = 0; n < signal.length; n++) {
				time[n] = n*dt;
			}

			// Close the WAV file
			wavFile.close();
		}
		catch (Exception e) {
			System.err.println(e);
		}
		
		//////////////////////////////////////		
		// 2. Padding
		//////////////////////////////////////

		// Calcuate the next highest power of 2 after the signal length
		double[] paddedSignal = new double[(int)Math.pow(2,((int)log(signal.length,2)+1))];
			
		// Loop to copy and pad the signal
		for (int i = 0; i < signal.length; i++) {
			paddedSignal [i] = signal [i] ;
		}
        for (int i = signal.length; i < paddedSignal.length; i++) {
            paddedSignal[i] = 0. ;
        }

		//////////////////////////////////////		
		// 3. FFT 
		//////////////////////////////////////

		// Create a FFT object and forward transform the signal
		FastFourierTransformer FFT = new FastFourierTransformer(DftNormalization.STANDARD);
		Complex[] signalFT = FFT.transform(paddedSignal,TransformType.FORWARD);

		//////////////////////////////////////		
		// 4. Power spectrum 
		//////////////////////////////////////
			
		// Create arrays for frequency and power
		double[] frequency = new double[signalFT.length/2];
		double[] power = new double[signalFT.length/2];
        double[] powerDistorted = new double[signalFT.length/2];
		
        // Calculate frequency and power and fill the arrays
		for (int m = 0; m < frequency.length; m++) {
			frequency[m] = m/dt/signalFT.length;
			power[m] = Math.pow(signalFT[m].getReal(),2)+Math.pow(signalFT[m].getImaginary(),2);
		}


		// Find maximum power
		double PMax = 0.0;
		for (int j = 0; j < power.length; j++) {
			if (power[j] > PMax) {
				PMax = power[j];
			}
		} 
        
        //////////////////////////////////////
        // 5. Finding the Fundamental Frequencies
        //////////////////////////////////////
        
        // detecting the number of fundamental frequencies to know how big fundFreq should be
        int numFundFreq = 0;
        
        for (int i = 0; i < power.length; i++) {
			if (power[i] >= arrayAvg(power)) {              
				numFundFreq++;
			}
		} 
		

        
        double[] fundFreq = new double [numFundFreq];   //stores the values of fundamental frequencies
        int[] fundFreqTracker = new int [numFundFreq];  //keeps track of the index of each fundamental frequency
    
        // finding the values of all the fundamental frequencies
        int count = 0;
        for (int i = 0; i < power.length; i++) {
            if (power[i] >= arrayAvg(power)) {              
				fundFreq[count] = frequency[i];
                fundFreqTracker[count] = i;     //this array keeps track of the index of the fundamental frequencies
                count++;
			}
		} 
        
        
        //////////////////////////////////////
        // 6. Adding The Harmonics / Distorting the Sound
        //////////////////////////////////////
        
        for (int i = 0; i < fundFreq.length; i++) {   //
            for (int j = 2; j < 8; j++) {         //this adds harmonics from the second harmonic to the 8th if that is within the range of our frequency spectrum   
                if (fundFreqTracker[i]*j < power.length) {
                    
                    signalFT[fundFreqTracker[i]*j] = signalFT[fundFreqTracker[i]].multiply(Math.sqrt(1./j));  //setting the real and imaginary values of the harmonic to a fraction of the fundamental
                    
                    signalFT[signalFT.length-1-fundFreqTracker[i]*j] = signalFT[fundFreqTracker[i]].multiply(Math.sqrt(1./j)); //doing the same operation on the corresponding index of the mirror image of the fourier transform to make sure that the symmetric quality of the fourier transform is conserved
                                                                                                                       }
                
            }
        }
        
        //////////////////////////////////////
        // 7. Inverse FFT
        //////////////////////////////////////                                                       
        // Putting the values of the transform into a new signal (which was padded)             
        Complex[] paddedSignal2 = FFT.transform(signalFT,TransformType.INVERSE);
		

        // Unpadding the signal, putting it into new array SignalDist
        for (int i = 0 ; i < signal.length ; i++) {
            signalDist[i] = paddedSignal2[i].getReal();
        }
        
        
        //converting the double array into an int array for writing
                                                                                                 
        int [] signalDistInt = new int [signalDist.length];
        
        for (int i = 0; i < signal.length ; i++) {
            signalDistInt[i] = (int)signalDist[i];
        }
        

        //////////////////////////////////////
        // 7. Writing the New Wav File
        //////////////////////////////////////
        
        try {		//parts taken from http://www.labbookpages.co.uk/audio/javaWavFiles.html
			
            double duration = 5.0;
            
            WavFile distortedSine = WavFile.newWavFile(new File(args[0]), 1, numFrames, 16, sampleRate);
			
            distortedSine.writeFrames(signalDist, numFrames);
            
			distortedSine.close();
		}
		catch (Exception e) {
			System.err.println(e);
		}
             

	}
        public static double log(double x, double base) {   //method that computes a log
            return (Math.log(x) / Math.log(base));
        }
        public static double arrayAvg(double [] array) {    //method that takes the average of an array
            double sum = 0;
            for (int i = 0; i < array.length; i++) {
                sum += array[i];
            }
            return (sum/array.length);
        }
        
}