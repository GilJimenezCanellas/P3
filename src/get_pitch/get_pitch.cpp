/// @file

#include <iostream>
#include <fstream>
#include <string.h>
#include <errno.h>
#include <algorithm> // for std::min and std::max

#include "wavfile_mono.h"
#include "pitch_analyzer.h"

#include "docopt.h"

#define FRAME_LEN   0.030 /* 30 ms. */
#define FRAME_SHIFT 0.015 /* 15 ms. */

using namespace std;
using namespace upc;

static const char USAGE[] = R"(
get_pitch - Pitch Estimator 

Usage:
    get_pitch [options] <input-wav> <output-txt>
    get_pitch (-h | --help)
    get_pitch --version

Options:
    -h, --help  Show this screen
    --version   Show the version of the project
    -t, --th_rlag=<th_rlag>   Set threshold for pitch estimation [default: 0.5]
    -r, --th_r1=<th_r1>   Set threshold for pitch estimation [default: 0.9]
    -z, --th_z=<th_z>   Threshold for zcr estimation [default: 30]
    -m, --th_m=<th_m>   Threshold for median filter [default: 3]
    -c, --th_c=<th_c>   Threshold for centre clipping [default: 0.0045]

Arguments:
    input-wav   Wave file with the audio signal
    output-txt  Output file: ASCII file with the result of the estimation:
                    - One line per frame with the estimated f0
                    - If considered unvoiced, f0 must be set to f0 = 0
)";

int main(int argc, const char *argv[]) {
	/// \TODO 
	///  Modify the program syntax and the call to *docopt()* in order to
	///  add options and arguments to the program.

  /// \DONE Program modified
    std::map<std::string, docopt::value> args = docopt::docopt(USAGE,
        {argv + 1, argv + argc},	// array of arguments, without the program name
        true,    // show help if requested
        "2.0");  // version string

	std::string input_wav = args["<input-wav>"].asString();
	std::string output_txt = args["<output-txt>"].asString();
  double th_rlag = std::stod(args["--th_rlag"].asString());
  double th_r1 = std::stod(args["--th_r1"].asString());
  double th_z = std::stod(args["--th_z"].asString());
  double th_m = std::stod(args["--th_m"].asString());
  double th_c = std::stod(args["--th_c"].asString());

  printf("TH_RLAG: %f\n", th_rlag);

  // Read input sound file
  unsigned int rate;
  vector<float> x;
  if (readwav_mono(input_wav, rate, x) != 0) {
    cerr << "Error reading input file " << input_wav << " (" << strerror(errno) << ")\n";
    return -2;
  }

  int n_len = rate * FRAME_LEN;
  int n_shift = rate * FRAME_SHIFT;

  // Define analyzer
  PitchAnalyzer analyzer(n_len, rate, PitchAnalyzer::RECT, 50, 500, th_rlag, th_r1, th_z);

  /// \TODO
  /// Preprocess the input signal in order to ease pitch estimation. For instance,
  /// central-clipping or low pass filtering may be used.

  /// \DONE center-clipping is done
  
  // Iterate for each frame and save values in f0 vector
  vector<float>::iterator iX;
  vector<float> f0;

  float alpha = 0.0047;
  for (iX = x.begin(); iX  < x.end(); iX++ ) {
    if (*iX < alpha && *iX > -alpha){ 
      *iX = 0;
    }
  }

  for (iX = x.begin(); iX + n_len < x.end(); iX = iX + n_shift) {
    float f = analyzer(iX, iX + n_len);
    f0.push_back(f);
  }

//******************
  // for (iX = x.begin(); iX + n_len < x.end(); iX = iX + n_shift) {
  //   // Calculate the maximum absolute value within the current frame
  //   float max_sample = *max_element(iX, iX + n_len, [](float a, float b) {
  //       return abs(a) < abs(b);
  //   });
  //   // Calculate the clip_threshold as half of the maximum absolute value
  //   float clip_threshold = 0.06 * abs(max_sample);

  //   // Apply central clipping
  //   for (int i = 0; i < n_len; ++i) {
  //     // Clip samples symmetrically around the center
  //     float sample = *(iX + i);
  //     if (abs(sample) < clip_threshold){
  //         *(iX + i) = 0;
  //     }
  //   }
    
  //   float f = analyzer(iX, iX + n_len);
  //   f0.push_back(f);
  // }
 // ****************
  /// \TODO
  /// Postprocess the estimation in order to supress errors. For instance, a median filter
  /// or time-warping may be used.

  /// \DONE median filter is done
//************
  // Apply median filter to the pitch estimation
  vector<float> filtered_f0;
  int filter_window_size = 3;

  for (int i = 0; i < f0.size(); ++i) {
      // Determine the range of indices for the current window
      int start_index = max(0, i - filter_window_size / 2);
      int end_index = min(static_cast<int>(f0.size()) - 1, i + filter_window_size / 2);
      
      // Create a copy of the pitch values within the window
      vector<float> window_pitch(f0.begin() + start_index, f0.begin() + end_index + 1);
      
      // Sort the window_pitch vector to find the median
      std::nth_element(window_pitch.begin(), window_pitch.begin() + window_pitch.size() / 2, window_pitch.end());
      
      // Get the median value
      float median_pitch = window_pitch[window_pitch.size() / 2];

      // This condition improves the pitch estimation, reduces the MSE fine errors,
      // but reduces the overall result
      // if (median_pitch - f0[i] < threshold1) {
      //   median_pitch = f0[i];
      // }
      
      // Store the median value in the filtered_f0 vector
      filtered_f0.push_back(median_pitch);
  }
//***********************

  // float aux = 0;
  // int k_wind = 1; ///window size
  // for(int i = 0; i < (int)f0.size(); i++){
  //   for(int j = 0; j < k_wind; ++j){
  //     if(i <= ((int)f0.size()-k_wind)){
  //       aux += f0[i + j]; 
  //     }else{
  //       aux = f0[i];
  //     }
  //   }
  //   f0[i] = aux/k_wind;
  //   aux=0;
  // }

  // Write f0 contour into the output file
  ofstream os(output_txt);
  if (!os.good()) {
    cerr << "Error reading output file " << output_txt << " (" << strerror(errno) << ")\n";
    return -3;
  }

  os << 0 << '\n'; //pitch at t=0
  for (iX = filtered_f0.begin(); iX != filtered_f0.end(); ++iX) 
    os << *iX << '\n';
  os << 0 << '\n';//pitch at t=Dur

  return 0;
}