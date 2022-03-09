//
// Fast Cross-Correlation algorithm
//
// Author: Francois Grondin
// Email: francois.grondin2@usherbrooke.ca
//
// Refer to the following paper for details:
//
// Grondin, F., Maheux, M.-A., Lauzon, J.-S. and Michaud, F., "Fast Cross-Correlation
// for TDoA Estimation on Small Aperture Microphone Arrays", arXiV
//

#ifndef __FCCPHAT_SIGNAL
#define __FCCPHAT_SIGNAL

    #ifdef __cplusplus
    extern "C" {
    #endif

    #include <stdlib.h>
    #include <string.h>
    #include <stdio.h>

    //
    // Signal that contains chunk of samples to fill frames
    //
    typedef struct hops_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per channel
        unsigned int hop_size;
        // Pointer to an array of arrays of samples
        float ** samples;

    } hops_obj;

    //
    // Signal that contains the frames in the frequency domain
    //
    typedef struct freqs_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per channel
        unsigned int frame_size;
        // Pointer to an array of arrays (size: channels_count) of frequency bins (size: frame_size / 2 + 1)
        float ** samples;

    } freqs_obj;

    //
    // Signal that contains the frames of the spatial covariance matrices
    //
    typedef struct covs_obj {

        // Number of channels
        unsigned int channels_count;
        // Number of samples per channel
        unsigned int frame_size;
        // Pointer to an array of arrays (size: channel_count * (channel_count-1) / 2) of cross-spectra (size: frame_size / 2 + 1)
        float ** samples;

    } covs_obj;

    //
    // Signal that contains the TDoAs with maximum correlation values
    //
    typedef struct corrs_obj {

        // Number of channels
        unsigned int channels_count;

        // The TDoA value before the TDoA with maximum correlation value (size: )
        float * taus_prev;
        // The correlation value that matches this TDoA
        float * ys_prev;

        // The TDoA value with maximum correlation value
        float * taus_max;
        // The correlation value that matches this TDoA
        float * ys_max;

        // The TDoA value after the TDoA with maximum correlation value
        float * taus_next;
        // The correlation value that matches this TDoA
        float * ys_next;

    } corrs_obj;

    //
    // Signal that contains the estimated TDoA values
    //
    typedef struct taus_obj {

        // Number of channels
        unsigned int channels_count;

        // Estimated TDoA values
        float * taus;
        // Corresponding maximum correlation amplitude
        float * ys;

    } taus_obj;

    //
    // Construct the hops object
    //
    // channels_count           Number of channels
    // hop_size                 Number of samples per channel
    //
    // (return)                 Pointer to the hops object
    //
    hops_obj * hops_construct(const unsigned int channels_count, const unsigned int hop_size);

    //
    // Destroy the hops object
    //
    // obj                      Pointer to the hops object
    //
    // (return)                 None
    //
    void hops_destroy(hops_obj * obj);

    //
    // Print the content of the hops object
    //
    // obj                      Pointer to the hops object
    //
    // (return)                 None, but prints content in console
    //
    void hops_printf(const hops_obj * obj);

    //
    // Construct the freqs object
    //
    // channels_count           Number of channels
    // frame_size               Number of samples per channel in each frame in the time-domain
    //
    // (return)                 Pointer to the freqs object
    //
    freqs_obj * freqs_construct(const unsigned int channels_count, const unsigned int frame_size);

    //
    // Destroy the freqs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None
    //
    void freqs_destroy(freqs_obj * obj);

    //
    // Print the content of the freqs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None, but prints content in console
    //
    void freqs_printf(const freqs_obj * obj);

    //
    // Construct the covs object
    //
    // channels_count           Number of channels
    // frame_size               Number of samples per channel in each frame in the time-domain
    //
    // (return)                 Pointer to the covs object
    //
    covs_obj * covs_construct(const unsigned int channels_count, const unsigned int frame_size);

    //
    // Destroy the covs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None
    //
    void covs_destroy(covs_obj * obj);

    //
    // Print the content of the covs object
    //
    // obj                      Pointer to the covs object
    //
    // (return)                 None, but prints content in console
    //
    void covs_printf(const covs_obj * obj);

    //
    // Construct the corrs object
    //
    // channels_count           Number of channels
    //
    // (return)                 Pointer to the corrs object
    //
    corrs_obj * corrs_construct(const unsigned int channels_count);

    //
    // Destroy the corrs object
    //
    // obj                      Pointer to the freqs object
    //
    // (return)                 None
    //
    void corrs_destroy(corrs_obj * obj);

    //
    // Print the content of the corrs object
    //
    // obj                      Pointer to the corrs object
    //
    // (return)                 None, but prints content in console
    //
    void corrs_printf(const corrs_obj * obj);

    //
    // Construct the taus object
    //
    // channels_count           Number of channels
    //
    // (return)                 Pointer to the taus object
    //
    taus_obj * taus_construct(const unsigned int channels_count);

    //
    // Destroy the taus object
    //
    // obj                      Pointer to the taus object
    //
    // (return)                 None
    //
    void taus_destroy(taus_obj * obj);

    //
    // Print the content of the taus object
    //
    // obj                      Pointer to the taus object
    //
    // (return)                 None, but prints content in console
    //
    void taus_printf(const taus_obj * obj);


    #ifdef __cplusplus
    }
    #endif

#endif
