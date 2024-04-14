/// @file

#include <iostream>
#include <math.h>
#include "pitch_analyzer.h"
#include <ffft/FFTReal.h>

using namespace std;

/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \TODO Compute the autocorrelation r[l]
      /// \DONE Autocorrelation computed
      /// \f[
      ///   r[l] = \frac{1}{N} \sum_n x[n][n+l]7
      /// \f]
      r[l] = 0.0f;
      for (unsigned int n = 0; n < x.size()-l; ++n) {
        r[l] += x[n] * x[n+l];
      }
      r[l] /= x.size();
    }

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }

  vector<float> PitchAnalyzer::cepstral_analysis(const vector<float> &x) const {

    // Perform the fft on the input signal
    ffft::FFTReal<float> fft(x.size());
    vector<float> x_fft(x.size()); // Complex FFT output
    fft.do_fft(&x_fft[0], &x[0]);
    
    // Compute the logarithm of the autocorrelation (cepstrum)
    for (unsigned int i = 0; i < x.size(); ++i) {
      float mod = sqrt(x_fft[i] * x_fft[i]);
      x_fft[i] = log(mod + 1e-10); // Add a small value to avoid log(0)
    }

    // Perform inverse Fourier transform to obtain the cepstrum
    vector<float> cepstrum(x.size());
    fft.do_ifft(&x_fft[0], &cepstrum[0]);
    fft.rescale(&cepstrum[0]);

    // Return the cepstrum
    return cepstrum;
  }

  pair<float, unsigned> PitchAnalyzer::get_results(const vector<float> &cepstrum) const {
    // Find the peak in the cepstrum (excluding the DC component)
    float max_val = -numeric_limits<float>::infinity();
    unsigned int max_idx = 0;
    for (unsigned int i = 10; i < cepstrum.size()/2; ++i) {
      if (cepstrum[i] > max_val) {
        max_val = cepstrum[i];
        max_idx = i;
      }
    }

    return make_pair(max_val, max_idx);
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    window.resize(frameLen);

    switch (win_type) {
    case HAMMING:
      /// \TODO Implement the Hamming window
      /// \DONE Hamming window implemented
      for (unsigned int n = 0; n < frameLen; ++n) {
        window[n] = 0.54f - 0.46f * cos(2.0f * M_PI * n / (frameLen - 1));
      }
      break;
    case RECT:
    default:
      window.assign(frameLen, 1);
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float max) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    if (max > 0.3) return false;
    return true;
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];

    // Perform cepstral analysis
    std::vector<float> cepstrum = cepstral_analysis(x);

    // Estimate pitch from cepstrum
    pair<float,unsigned> results = get_results(cepstrum);

    //******************** Here All the autocorrelation procedure starts
    // vector<float> r(npitch_max);

    // //Compute correlation
    // autocorrelation(x, r);

    // vector<float>::const_iterator iR = r.begin(), iRMax = iR;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.
    // iRMax = r.begin() + npitch_min;
    // for (iR = r.begin() + npitch_min; iR < r.begin() + npitch_max; iR++) {
    //   if (*iR > *iRMax) iRMax = iR;
    // }
    // unsigned int lag = iRMax - r.begin();

    // float pot = 10 * log10(r[0]);

    //********************** Here finishes

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif
    
    if (unvoiced(results.first))
      return 0;
    else
      return (float) samplingFreq/(float) results.second;
  }
}
