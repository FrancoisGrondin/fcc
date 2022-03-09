//
// Fast Cross-Correlation algorithm
//
// Authors: Francois Grondin, Marc-Antoine Maheux
// Emails: francois.grondin2@usherbrooke.ca, marc-antoine.maheux@usherbrooke.ca
//
// Refer to the following paper for details:
//
// Grondin, F., Maheux, M.-A., Lauzon, J.-S. and Michaud, F., "Fast Cross-Correlation
// for TDoA Estimation on Small Aperture Microphone Arrays", arXiV
//

#ifndef __FCCPHAT_SYSTEM
#define __FCCPHAT_SYSTEM

    #include "stdint.h"

    #ifdef __cplusplus
    extern "C" {
    #endif

    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <math.h>
    #include <fftw3.h>
    #include "signal.h"
    #include "const.h"

    //
    // Wave file header
    //
    // This structure holds the first 44 bytes used as a header
    // for the WAVE format.
    //
    typedef struct wav_header {

        // Chunk descriptor (corresponds to characters "RIFF")
        int8_t chunk_id[4];
        // Size of the file in bytes - 8
        int32_t chunk_size;
        // Contains the characters "WAVE"
        int8_t format[4];
        // Contains the characters "fmt "
        int8_t subchunk1_id[4];
        // Corresponds to 16 for PCM
        int32_t subchunk1_size;
        // PCM = 1 for linear quantization
        int16_t audio_format;
        // Number of channels
        int16_t num_channels;
        // Sample rate in samples/sec
        int32_t sample_rate;
        // Equals sample_rate * num_channels * bits_per_sample / 8
        int32_t byte_rate;
        // Equals num_channels * bits_per_sample / 8
        int16_t block_align;
        // Number of bits per sample
        int16_t bits_per_sample;
        // Contains the characters "data"
        int8_t subchunk2_id[4];
        // Equals num_samples * num_channels * bits_per_sample / 8
        int32_t subchunk2_size;

    } __attribute__((packed)) wav_header;

    //
    // Wave object that handles reading a WAVE file
    //
    typedef struct wav_obj {

        // Sample rate in samples/sec
        unsigned int sample_rate;
        // Number of channels
        unsigned int num_channels;
        // Number of bits per sample
        unsigned int bits_per_sample;
        // Number of samples read per channel each time
        unsigned int hop_size;

        // File pointer to read from the wave file
        FILE * file_pointer;
        // Buffer to load the samples (here we assume 16 bits per sample)
        int16_t * buffer;

    } wav_obj;

    //
    // STFT object to perform Short-Time Fourier Transform
    //
    typedef struct stft_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per frame to perform FFT
        unsigned int frame_size;
        // Hop size between frames in samples
        unsigned int hop_size;

        // This pointer refers to an array with the window's coefficients
        // (size = frame_size)
        float * window;
        // This pointer refers to an array of frames, each one containing an array of samples
        // (size = num_channels, and for each frame, size = frame_size)
        float ** frames;

        // This is the FFTW plan
        fftwf_plan fft;
        // This array points to the time-domain samples to perform the FFT
        // (size = frame_size)
        float * frame_real;
        // This array points to the frequency-domain samples obtained after FFT
        // (size = frame_size/2+1)
        fftwf_complex * frame_complex;

    } stft_obj;

    //
    // Spatial Covariance Matrix object to estimate the cross-spectra and
    // perform phase transform to normalize the amplitude in the frequency domain
    //
    typedef struct scmphat_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples in the time-domain prior to the STFT
        unsigned int frame_size;
        // Smoothing parameter to average over time
        // (alpha is between [0,1], and alpha=1 means no smoothing at all)
        float alpha;
        // Method to compute the result: if set to 'c', it will normalize each channel
        // before computing the cross-spectrum, and if set to 'p' it will compute the
        // cross-spectrum and then normalize
        char method;

        // Array of frames in the frequency domain with the cross-spectra
        float ** cross_spectrum;

    } scmphat_obj;

    //
    // Generalized Cross-Correlation object to estimate the TDoA
    //
    typedef struct gcc_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples in the time-domain prior to the STFT
        unsigned int frame_size;
        // Maximum magnitude for the TDoA (will be between [-tau_max,+tau_max])
        unsigned int tau_max;
        // Interpolation rate for the IFFT (r=1, r=2, r=4, ...)
        unsigned int interpolation_rate;

        // This is the FFTW plan
        fftwf_plan ifft;
        // This array points to the frequency-domain samples obtained after FFT
        // (size = frame_size/2+1)
        fftwf_complex * frame_complex;
        // This array points to the time-domain samples to perform the FFT
        // (size = frame_size)
        float * frame_real;

        // Array to hold only values in the interval [-tau_max,+tau_max]
        float * cropped_values;

    } gcc_obj;

    //
    // Fast Cross-Correlation object to estimate the TDoA
    //
    typedef struct fcc_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples in the time-domain prior to the STFT
        unsigned int frame_size;
        // Maximum magnitude for the TDoA (will be between [-tau_max,+tau_max])
        unsigned int tau_max;

        // Number of bases
        unsigned int K;
        // Number of TDoA candidates
        unsigned int L;
        // Array of bases
        float ** bases;
        // Array of keys for each TDoA
        float ** dicts;

        // This array holds the folded values with addition
        float * x_add;
        // This array holds the folded values with subtraction
        float * x_sub;
        // This array holds the projected vector
        float * z;
        // This array holds the power for each TDoA candidate
        float * y;

    } fcc_obj;

    //
    // Quadratric interpolation object
    //
    typedef struct quadinterp_obj {

        // Number of channels
        unsigned int channels_count;

    } quadinterp_obj;

    //
    // CSV object to write results to file
    //
    typedef struct csv_obj {

        // File pointer
        FILE * file_pointer;

        // Number of channels
        unsigned int channels_count;
        // The frame index
        unsigned int frame_index;

    } csv_obj;

    //
    // Construct the wave object
    //
    // file_name                Path of the wave file
    // hop_size                 Number of samples per channel to read at once
    //
    // (return)                 Pointer to the wave object
    //
    wav_obj * wav_construct(const char * file_name, const unsigned int hop_size);

    //
    // Destroy the wave object
    //
    // obj                      Pointer to the wave object
    //
    // (return)                 None
    //
    void wav_destroy(wav_obj * obj);

    //
    // Read samples from the wave object
    //
    // obj                  Pointer to the wave object
    // hops                 Pointer to the hops object that receives samples
    //
    // (return)             Returns -1 if the end of file reached, 0 otherwise
    //
    int wav_read(wav_obj * obj, hops_obj * hops);

    //
    // Construct the stft object
    //
    // channels_count       Number of channels
    // frame_size           Number of samples per frame
    // hop_size             Number of samples between frame
    //
    // (return)             Pointer of the stft object
    //
    stft_obj * stft_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int hop_size);

    //
    // Destroy the stft object
    //
    // obj                  Pointer to the stft object
    //
    // (return)             None
    //
    void stft_destroy(stft_obj * obj);

    //
    // Compute the Short-Time Fourier Transform
    //
    // obj                  Pointer to the stft object
    // hops                 Pointer to the hops object with new samples
    // freqs                Pointer to the freqs object with the results
    //
    // (return)             Returns 0 if no error
    //
    int stft_call(stft_obj * obj, const hops_obj * hops, freqs_obj * freqs);

    //
    // Construct the scmphat object
    //
    // channels_count       Number of channels
    // frame_size           Number of samples per frame
    // alpha                Smoothing factor
    //
    // (return)             Pointer to the scmphat object
    //
    scmphat_obj * scmphat_construct(const unsigned int channels_count, const unsigned int frame_size, const float alpha);

    //
    // Destroy the scmphat object
    //
    // obj                  Pointer to the scmphat object
    //
    // (return)             None
    //
    void scmphat_destroy(scmphat_obj * obj);

    //
    // Compute the Spatial Covariance Matrix and normalize
    //
    // obj                  Pointer to the scmphat object
    // freqs                Pointer to the freqs object with stft samples
    // covs                 Pointer to the covs object with scm results
    //
    // (return)             Returns 0 if no error
    //
    int scmphat_call(scmphat_obj * obj, const freqs_obj * freqs, covs_obj * covs);

    //
    // Construct the gcc object
    //
    // channels_count       Number of channels
    // frame_size           Number of samples per frame
    // tau_max              Maximum magnitude of the TDoA
    // interpolation_rate   Interpolation rate (r=1, r=2, r=4, ...)
    //
    // (return)             Pointer to the gcc object
    //
    gcc_obj * gcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max, const unsigned int interpolation_rate);

    //
    // Destroy the gcc object
    //
    // obj                  Pointer to the gcc object
    //
    // (return)             None
    //
    void gcc_destroy(gcc_obj * obj);

    //
    // Compute Generalized Cross-Correlation
    //
    // obj                  Pointer to the gcc object
    // covs                 Pointer to the covs object that holds the SCM
    // corrs                Pointer to the corrs object with the cross-correlation in time-domain
    //
    // (return)             Returns 0 if no error
    //
    int gcc_call(gcc_obj * obj, const covs_obj * covs, corrs_obj * corrs);

    //
    // Construct the fcc object
    //
    // channels_count       Number of channels
    // frame_size           Number of samples per frame
    // tau_max              Maximum magnitude of the TDoA
    //
    // (return)             Pointer to the fcc object
    //
    fcc_obj * fcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max);

    //
    // Destroy the fcc object
    //
    // obj                  Pointer to the fcc object
    //
    // (return)             None
    //
    void fcc_destroy(fcc_obj * obj);

    //
    // Compute Fast Cross-Correlation
    //
    // obj                  Pointer to the fcc object
    // covs                 Pointer to the covs object that holds the SCM
    // corrs                Pointer to the corrs object with the cross-correlation in time-domain
    //
    // (return)             Returns 0 if no error
    //
    int fcc_call(fcc_obj * obj, const covs_obj * covs, corrs_obj * corrs);

    //
    // Construct the quadinterp object
    //
    // channels_count       Number of channels
    //
    // (return)             Pointer to the quadinterp object
    //
    quadinterp_obj * quadinterp_construct(const unsigned int channels_count);

    //
    // Destroy the quadinterp object
    //
    // obj                  Pointer to the quadinterp object
    //
    // (return)             None
    //
    void quadinterp_destroy(quadinterp_obj * obj);

    //
    // Compute the quadratic interpolation
    //
    // obj                  Pointer to the quadinterp object
    // corrs                Pointer to the corrs object with the cross-correlation in time-domain
    // taus                 Pointer to the taus object with the tdoas
    //
    // (return)             Returns 0 if no error
    //
    int quadinterp_call(quadinterp_obj * obj, const corrs_obj * corrs, taus_obj * taus);

    //
    // Construct the csv object
    //
    // file_name            Path of the file name to write the CSV content
    // channels_count       Number of channels
    //
    // (return)             Pointer to the csv object
    //
    csv_obj * csv_construct(const char * file_name, const unsigned int channels_count);

    //
    // Destroy the csv object
    //
    // obj                  Pointer to the csv object
    //
    // (return)             None
    //
    void csv_destroy(csv_obj * obj);

    //
    // Write the TDoAs to the CSV file
    //
    // obj                  Pointer to the csv object
    // taus                 Pointer to the taus object which contains the TDoAs
    //
    // (return)             Returns 0 if no error
    //
    int csv_write(csv_obj * obj, taus_obj * taus);

    #ifdef __cplusplus
    }
    #endif

#endif
