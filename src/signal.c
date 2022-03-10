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

#include <fccphat/signal.h>
#include <fccphat/simd.h>

//
// Construct the hops object
//
// channels_count           Number of channels
// hop_size                 Number of samples per channel
//
// (return)                 Pointer to the hops object
//
hops_obj * hops_construct(const unsigned int channels_count, const unsigned int hop_size) {

    hops_obj * obj;
    unsigned int channel_index;

    // Allocate memory for object
    obj = (hops_obj *) malloc(sizeof(hops_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->hop_size = hop_size;

    // Create an array, with a frame with hop_size samples for each channel
    obj->samples = (float **) malloc(sizeof(float *) * channels_count);
    for (channel_index = 0; channel_index < channels_count; channel_index++) {
        obj->samples[channel_index] = (float *) malloc(sizeof(float) * hop_size);
        memset(obj->samples[channel_index], 0x00, sizeof(float) * hop_size);
    }

    // Return pointer to object
    return obj;

}

//
// Destroy the hops object
//
// obj                      Pointer to the hops object
//
// (return)                 None
//
void hops_destroy(hops_obj * obj) {

    unsigned int channel_index;

    // Deallocate memory for frames
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        free((void *) obj->samples[channel_index]);
    }
    free((void *) obj->samples);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Print the content of the hops object
//
// obj                      Pointer to the hops object
//
// (return)                 None, but prints content in console
//
void hops_printf(const hops_obj * obj) {

    unsigned int channel_index;
    unsigned int sample_index;

    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

        for (sample_index = 0; sample_index < obj->hop_size; sample_index++){

            printf("[%02u]-(%03u): %1.3f\n", channel_index, sample_index, obj->samples[channel_index][sample_index]);

        }

    }

}

//
// Construct the freqs object
//
// channels_count           Number of channels
// frame_size               Number of samples per channel in each frame in the time-domain
//
// (return)                 Pointer to the freqs object
//
freqs_obj * freqs_construct(const unsigned int channels_count, const unsigned int frame_size) {

    freqs_obj * obj;
    unsigned int channel_index;

    // Allocate memory for object
    obj = (freqs_obj *) malloc(sizeof(freqs_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;

    // Create an array, with a frame of frame_size/2+1 complex numbers for each channel
    obj->samples = (float **) malloc(sizeof(float *) * channels_count);
    for (channel_index = 0; channel_index < channels_count; channel_index++) {
        obj->samples[channel_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
        memset(obj->samples[channel_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
    }

    // Return pointer to object
    return obj;

}

//
// Destroy the freqs object
//
// obj                      Pointer to the freqs object
//
// (return)                 None
//
void freqs_destroy(freqs_obj * obj) {

    unsigned int channel_index;

    // Deallocate memory for frames
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        free((void *) obj->samples[channel_index]);
    }
    free((void *) obj->samples);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Print the content of the freqs object
//
// obj                      Pointer to the freqs object
//
// (return)                 None, but prints content in console
//
void freqs_printf(const freqs_obj * obj) {

    unsigned int channel_index;
    unsigned int bin_index;

    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

        for (bin_index = 0; bin_index < obj->frame_size/2+1; bin_index++){

            printf("[%02u]-(%03u): (%+1.3f , %+1.3f)\n", channel_index, bin_index, obj->samples[channel_index][bin_index*2+0], obj->samples[channel_index][bin_index*2+1]);

        }

    }

}

//
// Construct the covs object
//
// channels_count           Number of channels
// frame_size               Number of samples per channel in each frame in the time-domain
//
// (return)                 Pointer to the covs object
//
covs_obj * covs_construct(const unsigned int channels_count, const unsigned int frame_size) {

    covs_obj * obj;
    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    // Allocate memory for object
    obj = (covs_obj *) malloc(sizeof(covs_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;

    // For each pair of microphones, create a frame with (frame_size/2+1) complex numbers
    pair_index = 0;
    obj->samples = (float **) malloc(sizeof(float *) * channels_count * (channels_count-1) / 2);
    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
        for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
#ifdef FCCPHAT_USE_SIMD
            obj->samples[pair_index] = (float *) aligned_alloc(FLOAT_SIMD_ALIGNMENT, sizeof(float) * (frame_size/2+1) * 2);
#else
            obj->samples[pair_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
#endif

            memset(obj->samples[pair_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
            pair_index++;
        }
    }

    // Return pointer to object
    return obj;

}

//
// Destroy the covs object
//
// obj                      Pointer to the freqs object
//
// (return)                 None
//
void covs_destroy(covs_obj * obj) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    // Deallocate memory for frames
    pair_index = 0;
    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
        for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
            free((void *) obj->samples[pair_index]);
            pair_index++;
        }
    }
    free((void *) obj->samples);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Print the content of the covs object
//
// obj                      Pointer to the covs object
//
// (return)                 None, but prints content in console
//
void covs_printf(const covs_obj * obj) {

    unsigned int pair_index;
    unsigned int bin_index;

    for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count-1)/2); pair_index++) {

        for (bin_index = 0; bin_index < obj->frame_size/2+1; bin_index++){

            printf("[%02u]-(%03u): (%+1.3f , %+1.3f)\n", pair_index, bin_index, obj->samples[pair_index][bin_index*2+0], obj->samples[pair_index][bin_index*2+1]);

        }

    }

}

//
// Construct the corrs object
//
// channels_count           Number of channels
//
// (return)                 Pointer to the corrs object
//
corrs_obj * corrs_construct(const unsigned int channels_count) {

    corrs_obj * obj;

    // Allocate memory for object
    obj = (corrs_obj *) malloc(sizeof(corrs_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;

    // Allocate memory for TDoAs and corresponding cross-correlation values for each pair of microphones
    obj->taus_prev = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->taus_prev, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
    obj->ys_prev = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->ys_prev, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
    obj->taus_max = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->taus_max, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
    obj->ys_max = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->ys_max, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
    obj->taus_next = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->taus_next, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
    obj->ys_next = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->ys_next, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);

    // Return pointer to object
    return obj;

}

//
// Destroy the corrs object
//
// obj                      Pointer to the freqs object
//
// (return)                 None
//
void corrs_destroy(corrs_obj * obj) {

    // Deallocate memory for TDoAs and corresponding cross-correlation values
    free((void *) obj->taus_prev);
    free((void *) obj->ys_prev);
    free((void *) obj->taus_max);
    free((void *) obj->ys_max);
    free((void *) obj->taus_next);
    free((void *) obj->ys_next);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Print the content of the corrs object
//
// obj                      Pointer to the corrs object
//
// (return)                 None, but prints content in console
//
void corrs_printf(const corrs_obj * obj) {

    unsigned int pair_index;

    for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count-1)/2); pair_index++) {

        printf("[%02u] taus_prev=%+1.3f ys_prev=%+1.3f taus_max=%+1.3f ys_max=%+1.3f taus_next=%+1.3f ys_next=%+1.3f\n",
            pair_index, obj->taus_prev[pair_index], obj->ys_prev[pair_index], obj->taus_max[pair_index], obj->ys_max[pair_index], obj->taus_next[pair_index], obj->ys_next[pair_index]);

    }

}

//
// Construct the taus object
//
// channels_count           Number of channels
//
// (return)                 Pointer to the taus object
//
taus_obj * taus_construct(const unsigned int channels_count) {

    taus_obj * obj;

    // Allocate memory for object
    obj = (taus_obj *) malloc(sizeof(taus_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;

    // Allocate memory for the estimated TDoAs and corresponding cross-correlation values
    obj->taus = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->taus, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
    obj->ys = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
    memset(obj->ys, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);

    // Return pointer to object
    return obj;

}

//
// Destroy the taus object
//
// obj                      Pointer to the taus object
//
// (return)                 None
//
void taus_destroy(taus_obj * obj) {

    // Deallocate memory for TDoAs and corresponding cross-correlation values
    free((void *) obj->taus);
    free((void *) obj->ys);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Print the content of the taus object
//
// obj                      Pointer to the taus object
//
// (return)                 None, but prints content in console
//
void taus_printf(const taus_obj * obj) {

    unsigned int pair_index;

    for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count - 1)/2); pair_index++) {

        printf("[%02u] taus=%1.3f ys=%1.3f\n", pair_index, obj->taus[pair_index], obj->ys[pair_index]);

    }

}
