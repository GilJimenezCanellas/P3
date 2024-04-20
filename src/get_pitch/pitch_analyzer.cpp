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

  void PitchAnalyzer::mdf(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \TODO Compute the autocorrelation r[l]
      /// \DONE MDF computed
      /// \f[
      ///   r[l] = \frac{1}{N} \sum_n x[n]-[n+l]7
      /// \f]
      r[l] = 0.0f;
      for (unsigned int n = 0; n < x.size()-l; ++n) {
        r[l] += abs(x[n] - x[n+l]);
      }
      r[l] /= x.size();
    }

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 
  }

  vector<float> PitchAnalyzer::cepstral_analysis(const vector<float> &x) const {

    vector<float> x_corrected = x;

    // Add zeros to have a length multiple of 2
    size_t next_power_of_2 = pow(2, ceil(log2(x_corrected.size())));
    size_t num_zeros = next_power_of_2 - x_corrected.size();
    x_corrected.insert(x_corrected.end(), num_zeros, 0.0f);

    // Perform the fft on the input signal
    ffft::FFTReal<float> fft(x_corrected.size());
    float x_[x_corrected.size()];
    for (unsigned int i = 0; i < x_corrected.size(); ++i) 
      x_[i] = x_corrected[i];
    float x_fft[x_corrected.size()]; // Complex FFT output
    fft.do_fft(x_fft, x_);
    
    // Compute the logarithm of the autocorrelation (cepstrum)
    for (unsigned int i = 0; i < x_corrected.size(); ++i) {
      float mod = sqrt(x_fft[i] * x_fft[i]);
      x_fft[i] = log(mod + 1e-10); // Add a small value to avoid log(0)
    }

    // Perform inverse Fourier transform to obtain the cepstrum
    float cepstrum_[x_corrected.size()];
    fft.do_ifft(x_fft, cepstrum_);
    fft.rescale(cepstrum_);

    vector<float> cepstrum;
    cepstrum.resize(x_corrected.size());
    for (int i = 0; i < x_corrected.size(); ++i)
      cepstrum[i] = cepstrum_[i];

    return cepstrum;
  }

  tuple<float, unsigned, float> PitchAnalyzer::get_results(const vector<float> &cepstrum) const {
    // Find the peak in the cepstrum (excluding the DC component)
    float max_val = -10000;
    unsigned int max_idx = 0;
    float max_val_zero = 0;
    for (unsigned int i= 0; i < 10; ++i) {
      if (cepstrum[i] > max_val_zero) {
        max_val_zero = cepstrum[i];
      }
    }
    for (unsigned int i = 40; i < 400; ++i) {
      if (cepstrum[i] > max_val) {
        max_val = cepstrum[i];
        max_idx = i;
      }
    }

    return make_tuple(max_val, max_idx, max_val_zero);
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

  bool PitchAnalyzer::unvoiced(int zcr, float r1norm, float rmaxnorm, float pot) const {
    /// \TODO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
    // if (rmaxnorm > 0.34 && r1norm > 0.9 && zcr < 150) return false;
    // int count = 0;

    // if (rmaxnorm > 0.35) count++;
    // if (r1norm > 0.87) count++;
    // if (pot > threshold1) count++;

    // if (count >= 2) return false;

    if (rmaxnorm > 0.35 && r1norm > 0.87) return false;
    return true;
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i)
      x[i] *= window[i];
    
    // Perform cepstral analysis
    // std::vector<float> cepstrum = cepstral_analysis(x);

    // Estimate pitch from cepstrum
    // tuple<float, unsigned, float> results = get_results(cepstrum);

    //******************** Here All the autocorrelation procedure starts
    vector<float> r(npitch_max);
    vector<float> d(npitch_max);

    //Compute correlation
    autocorrelation(x, r);
    // mdf(x, d);

    int zcr = 0;
    // for (int i = 1; i < x.size(); i++) {
    //     zcr += (x[i] * x[i-1] < 0);
    // }

    vector<float>::const_iterator iR = r.begin(), iRMax = iR;
    // vector<float>::const_iterator iD = d.begin(), iDMin = iD;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.
    iRMax = r.begin() + npitch_min;
    for (iR = r.begin() + npitch_min; iR < r.begin() + npitch_max; iR++) {
      if (*iR > *iRMax) iRMax = iR;
    }
    unsigned int lag = iRMax - r.begin();

    float pot = 10 * log10(r[0]);

    // iDMin = d.begin() + npitch_min;
    // for (iD = d.begin() + npitch_min; iD < d.begin() + npitch_max; iD++) {
    //   if (*iD < *iDMin) iDMin = iD;
    // }
    // unsigned int lag_d = iDMin - d.begin();

    //********************** Here finishes

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 0
    if (r[0] > 0.0F)
      cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif
    
    if (unvoiced(zcr, r[1]/r[0], r[lag]/r[0], pot))
      return 0;
    else
      return (float) samplingFreq/(float) lag;
  }
}
