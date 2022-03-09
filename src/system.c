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

#include <fccphat/system.h>

//
// Construct the wave object
//
// file_name                Path of the wave file
// hop_size                 Number of samples per channel to read at once
//
// (return)                 Pointer to the wave object
//
wav_obj * wav_construct(const char * file_name, const unsigned int hop_size) {

    wav_obj * obj;
    wav_header hdr;
    size_t rtn;

    // Allocate memory for object
    obj = (wav_obj *) malloc(sizeof(wav_obj));

    // Open the wave file to be read
    obj->file_pointer = fopen(file_name, "rb");
    if (obj->file_pointer == NULL) {
        printf("File %s does not exists\n", file_name);
        exit(-1);
    }

    // Load the first 44 bytes that are the header of the file
    rtn = fread(&hdr, sizeof(char), sizeof(wav_header), obj->file_pointer);
    if (rtn != sizeof(wav_header)) {
        printf("Invalid WAV file.\n");
        exit(-1);
    }

    // Copy relevant parameters
    obj->sample_rate = hdr.sample_rate;
    obj->num_channels = hdr.num_channels;
    obj->bits_per_sample = hdr.bits_per_sample;
    obj->hop_size = hop_size;

    // Here we only support signed 16-bit samples
    if (obj->bits_per_sample != 16) {
        printf("Unsupported PCM format.\n");
        exit(-1);
    }

    // Create a buffer of shorts with a total of (hop_size * num_channels) samples
    obj->buffer = (short *) malloc(sizeof(short) * obj->hop_size * obj->num_channels);

    // Return the pointer to the wave object
    return obj;

}

//
// Destroy the wave object
//
// obj                      Pointer to the wave object
//
// (return)                 None
//
void wav_destroy(wav_obj * obj) {

    // First close the wave file
    fclose(obj->file_pointer);
    // Then free the buffer that held all the samples
    free((void *) obj->buffer);

    // Finally deallocate memory for object
    free((void *) obj);

}

//
// Read samples from the wave object
//
// obj                  Pointer to the wave object
// hops                 Pointer to the hops object that receives samples
//
// (return)             Returns -1 if the end of file reached, 0 otherwise
//
int wav_read(wav_obj * obj, hops_obj * hops) {

    size_t rtn;
    unsigned int channel_index;
    unsigned int sample_index;

    // Read from the file a total of 2 * hop_size * num_channels bytes
    rtn = fread(obj->buffer, obj->bits_per_sample/8, obj->hop_size * obj->num_channels, obj->file_pointer);

    // Then load each sample in the buffer (samples are interleaved, and are signed 16-bit)
    // Then are divided by 32768 to be normalized between -1 and 1
    for (sample_index = 0; sample_index < obj->hop_size; sample_index++) {
        for (channel_index = 0; channel_index < obj->num_channels; channel_index++) {
            hops->samples[channel_index][sample_index] = ((float) obj->buffer[sample_index * obj->num_channels + channel_index]) / 32768.0;
        }
    }

    // If the number of elements read is consistent, return 0 (meaning there are still elements to read)
    // otherwise return -1, which means we've reached the end of the file
    return (rtn == obj->hop_size * obj->num_channels) ? 0:-1;

}

//
// Construct the stft object
//
// channels_count       Number of channels
// frame_size           Number of samples per frame
// hop_size             Number of samples between frame
//
// (return)             Pointer of the stft object
//
stft_obj * stft_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int hop_size) {

    stft_obj * obj;
    unsigned int channel_index;
    unsigned int sample_index;

    // Allocate memory for object
    obj = (stft_obj *) malloc(sizeof(stft_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;
    obj->hop_size = hop_size;

    // Create a Hann window
    obj->window = (float *) malloc(sizeof(float) * frame_size);
    for (sample_index = 0; sample_index < frame_size; sample_index++) {
        obj->window[sample_index] = powf(sinf((float) sample_index * M_PI / ((float) (frame_size - 1))), 2.0);
    }

    // Allocate memory for each frame. There is one frame of frame_size samples per channel.
    obj->frames = (float **) malloc(sizeof(float *) * channels_count);
    for (channel_index = 0; channel_index < channels_count; channel_index++) {
        obj->frames[channel_index] = (float *) malloc(sizeof(float) * frame_size);
        memset((void *) obj->frames[channel_index], 0x00, sizeof(float) * frame_size);
    }

    // Allocate memory to perform the FFT using FFTW
    obj->frame_real = (float *) fftwf_malloc(sizeof(float) * frame_size);
    obj->frame_complex = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * (frame_size/2+1));
    obj->fft = fftwf_plan_dft_r2c_1d(frame_size, obj->frame_real, obj->frame_complex, FFTW_EXHAUSTIVE);

    // Return the pointer to the stft object
    return obj;

}

//
// Destroy the stft object
//
// obj                  Pointer to the stft object
//
// (return)             None
//
void stft_destroy(stft_obj * obj) {

    unsigned int channel_index;

    // Deallocate memory for window
    free((void *) obj->window);

    // Deallocate memory for each frame, and then clear list of pointers
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
        free((void *) obj->frames[channel_index]);
    }
    free((void *) obj->frames);

    // Deallocate memory for FFTW
    fftwf_free(obj->frame_real);
    fftwf_free(obj->frame_complex);
    fftwf_destroy_plan(obj->fft);
    fftwf_cleanup();

    // Finally, deallocate memory for the object
    free((void *) obj);

}

//
// Compute the Short-Time Fourier Transform
//
// obj                  Pointer to the stft object
// hops                 Pointer to the hops object with new samples
// freqs                Pointer to the freqs object with the results
//
// (return)             Returns 0 if no error
//
int stft_call(stft_obj * obj, const hops_obj * hops, freqs_obj * freqs) {

    unsigned int channel_index;
    unsigned int sample_index;
    unsigned int bin_index;

    // Load the samples for each channel
    for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

        //
        // Shift samples to the left:
        //
        // The initial array
        //
        // |AAAAAAAAAAAAAA|BBBBBBBBBBBBBBBBBBBBBBBBBBB|
        //  <--hop_size--> <--(frame_size-hop_size)-->
        //
        // becomes
        //
        // |BBBBBBBBBBBBBBBBBBBBBBBBBBB|00000000000000|
        //  <--(frame_size-hop_size)--> <--hop_size-->
        //
        memmove(&(obj->frames[channel_index][0]), &(obj->frames[channel_index][obj->hop_size]), sizeof(float) * (obj->frame_size - obj->hop_size));

        //
        // Then fill the right side with new samples
        //
        // So the array
        //
        // |BBBBBBBBBBBBBBBBBBBBBBBBBBB|00000000000000|
        //  <--(frame_size-hop_size)--> <--hop_size-->
        //
        // becomes
        //
        // |BBBBBBBBBBBBBBBBBBBBBBBBBBB|CCCCCCCCCCCCCC|
        //  <--(frame_size-hop_size)--> <--hop_size-->
        //
        memcpy(&(obj->frames[channel_index][obj->frame_size - obj->hop_size]), &(hops->samples[channel_index][0]), sizeof(float) * obj->hop_size);

        // Then loop for each sample in the frame and multiply it by the window, and save result to array to perform FFT
        for (sample_index = 0; sample_index < obj->frame_size; sample_index++) {
            obj->frame_real[sample_index] = obj->window[sample_index] * obj->frames[channel_index][sample_index];
        }

        // Perform FFT (result goes in frame_complex)
        fftwf_execute(obj->fft);

        // Extract the results and copy them to the freqs object
        for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {
            freqs->samples[channel_index][bin_index * 2 + 0] = obj->frame_complex[bin_index][0];
            freqs->samples[channel_index][bin_index * 2 + 1] = obj->frame_complex[bin_index][1];
        }

    }

    // By default return 0 when returns without error
    return 0;

}

//
// Construct the scmphat object
//
// channels_count       Number of channels
// frame_size           Number of samples per frame
// alpha                Smoothing factor
//
// (return)             Pointer to the scmphat object
//
scmphat_obj * scmphat_construct(const unsigned int channels_count, const unsigned int frame_size, const float alpha) {

    scmphat_obj * obj;
    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;
    unsigned int channel_index;

    // Allocate memory for the object
    obj = (scmphat_obj *) malloc(sizeof(scmphat_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;
    obj->alpha = alpha;

    // If we take only the instantaneous cross-correlation (alpha = 1)
    // and there are more than 2 channels, it is faster to perform PHAT
    // normalization on each channel and then do the pairwise multiplication
    obj->method = 'p';
    if (alpha == 1.0f) {
        if (channels_count > 2) {
            obj->method = 'c';
        }
    }

    // If pairwise approach, then assign the cross spectrum arrays
    if (obj->method == 'p') {

        pair_index = 0;
        obj->cross_spectrum = (float **) malloc(sizeof(float *) * (channels_count * (channels_count-1) / 2));
        for (channel_index1 = 0; channel_index1 < channels_count; channel_index1++) {
            for (channel_index2 = (channel_index1+1); channel_index2 < channels_count; channel_index2++) {
                obj->cross_spectrum[pair_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
                memset((void *) obj->cross_spectrum[pair_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
                pair_index++;
            }
        }

    }
    // If channelwise approach, then assign the cross spectrum arrays
    if (obj->method == 'c') {

        obj->cross_spectrum = (float **) malloc(sizeof(float *) * channels_count);
        for (channel_index = 0; channel_index < channels_count; channel_index++) {
            obj->cross_spectrum[channel_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
            memset((void *) obj->cross_spectrum[channel_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
        }

    }

    // Return pointer to the object
    return obj;

}

//
// Destroy the scmphat object
//
// obj                  Pointer to the scmphat object
//
// (return)             None
//
void scmphat_destroy(scmphat_obj * obj) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;
    unsigned int channel_index;

    // If pairwise approach, deallocate memory accordingly
    if (obj->method == 'p') {

        pair_index = 0;
        for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
            for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
                free((void *) obj->cross_spectrum[pair_index]);
                pair_index++;
            }
        }
        free((void *) obj->cross_spectrum);

    }
    // Or deallocate in the channelwise approach
    if (obj->method == 'c') {

        for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
            free((void *) obj->cross_spectrum[channel_index]);
        }
        free((void *) obj->cross_spectrum);

    }

    // Deallocate memory for object
    free((void *) obj);

}

//
// Compute the Spatial Covariance Matrix and normalize
//
// obj                  Pointer to the scmphat object
// freqs                Pointer to the freqs object with stft samples
// covs                 Pointer to the covs object with scm results
//
// (return)             Returns 0 if no error
//
int scmphat_call(scmphat_obj * obj, const freqs_obj * freqs, covs_obj * covs) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int bin_index;
    unsigned int pair_index;
    unsigned int channel_index;

    float dest_real;
    float dest_imag;
    float dest_magn;
    float src1_real;
    float src1_imag;
    float src2_real;
    float src2_imag;
    float src_real;
    float src_imag;
    float src_magn;

    // If pairwise approach
    if (obj->method == 'p') {

        // Compute the following for each pair of microphones i and j:
        //
        // 1) R_(i,j)(t,f) = (1 - alpha) * R_(i,j)(t-1,f) + alpha * X_i(t,f) * X_j(t,f)^*
        // 2) Rhat_(i,j)(t,f) = R_(i,j)(t,f) / |R_(i,j)(t,f)|

        pair_index = 0;

        for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

            for (channel_index2 = (channel_index1 + 1); channel_index2 < obj->channels_count; channel_index2++) {

                for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                    dest_real = obj->cross_spectrum[pair_index][bin_index*2+0];
                    dest_imag = obj->cross_spectrum[pair_index][bin_index*2+1];

                    src1_real = freqs->samples[channel_index1][bin_index*2+0];
                    src1_imag = freqs->samples[channel_index1][bin_index*2+1];
                    src2_real = freqs->samples[channel_index2][bin_index*2+0];
                    src2_imag = freqs->samples[channel_index2][bin_index*2+1];

                    dest_real *= (1.0 - obj->alpha);
                    dest_imag *= (1.0 - obj->alpha);

                    dest_real += obj->alpha * (src1_real * src2_real + src1_imag * src2_imag);
                    dest_imag += obj->alpha * (src1_imag * src2_real - src1_real * src2_imag);

                    obj->cross_spectrum[pair_index][bin_index*2+0] = dest_real;
                    obj->cross_spectrum[pair_index][bin_index*2+1] = dest_imag;

                    dest_magn = sqrtf(dest_real * dest_real + dest_imag * dest_imag);

                    dest_real /= (dest_magn + 1e-10);
                    dest_imag /= (dest_magn + 1e-10);

                    covs->samples[pair_index][bin_index*2+0] = dest_real;
                    covs->samples[pair_index][bin_index*2+1] = dest_imag;

                }

                pair_index++;

            }

        }

    }
    // If channelwise approach
    if (obj->method == 'c') {

        // Compute the following
        //
        // Rhat_(i,j)(t,f) = (X_i(t,f)/|X_i(t,f)|) * (X_j(t,f)^*/|X_j(t,f)|)

        for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

            for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                src_real = freqs->samples[channel_index][bin_index*2+0];
                src_imag = freqs->samples[channel_index][bin_index*2+1];
                src_magn = sqrtf(src_real * src_real + src_imag * src_imag);

                src_real /= (src_magn + 1e-10);
                src_imag /= (src_imag + 1e-10);

                obj->cross_spectrum[channel_index][bin_index*2+0] = src_real;
                obj->cross_spectrum[channel_index][bin_index*2+1] = src_imag;

            }

        }

        pair_index = 0;

        for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

            for (channel_index2 = (channel_index1 + 1); channel_index2 < obj->channels_count; channel_index2++) {

                for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                    src1_real = obj->cross_spectrum[channel_index1][bin_index*2+0];
                    src1_imag = obj->cross_spectrum[channel_index1][bin_index*2+1];
                    src2_real = obj->cross_spectrum[channel_index2][bin_index*2+0];
                    src2_imag = obj->cross_spectrum[channel_index2][bin_index*2+1];

                    dest_real = src1_real * src2_real + src1_imag * src2_imag;
                    dest_imag = src1_imag * src2_real - src1_real * src2_imag;

                    covs->samples[pair_index][bin_index*2+0] = dest_real;
                    covs->samples[pair_index][bin_index*2+1] = dest_imag;

                }

                pair_index++;

            }

        }

    }

    // By default return 0 when returns without error
    return 0;

}

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
gcc_obj * gcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max, const unsigned int interpolation_rate) {

    gcc_obj * obj;

    // Allocate memory for object
    obj = (gcc_obj *) malloc(sizeof(gcc_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;
    obj->tau_max = tau_max;
    obj->interpolation_rate = interpolation_rate;

    // Allocate memory for FFTW
    obj->frame_complex = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * (frame_size/2+1) * interpolation_rate);
    obj->frame_real = (float *) fftwf_malloc(sizeof(float) * frame_size * interpolation_rate);
    obj->ifft = fftwf_plan_dft_c2r_1d(frame_size * interpolation_rate, obj->frame_complex, obj->frame_real, FFTW_EXHAUSTIVE);

    // Create frame to hold only relevant interval of TDoAs
    obj->cropped_values = (float *) malloc(sizeof(float) * (tau_max * 2 * interpolation_rate + 1));

    // Return pointer to the object
    return obj;

}

//
// Destroy the gcc object
//
// obj                  Pointer to the gcc object
//
// (return)             None
//
void gcc_destroy(gcc_obj * obj) {

    // Deallocate memory for FFTW
    free((void *) obj->frame_complex);
    free((void *) obj->frame_real);
    fftwf_destroy_plan(obj->ifft);
    fftwf_cleanup();

    // Deallocate memory for array that held TDoAs
    free((void *) obj->cropped_values);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Compute Generalized Cross-Correlation
//
// obj                  Pointer to the gcc object
// covs                 Pointer to the covs object that holds the SCM
// corrs                Pointer to the corrs object with the cross-correlation in time-domain
//
// (return)             Returns 0 if no error
//
int gcc_call(gcc_obj * obj, const covs_obj * covs, corrs_obj * corrs) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;
    unsigned int bin_index;
    unsigned int tau_index;
    unsigned int max_index;

    pair_index = 0;

    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

        for (channel_index2 = (channel_index1 + 1); channel_index2 < obj->channels_count; channel_index2++) {

            // For each pair of microphones, compute the IFFT

            memset(obj->frame_complex, 0x00, sizeof(fftwf_complex) * (obj->frame_size/2+1) * obj->interpolation_rate);

            for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

                obj->frame_complex[bin_index][0] = covs->samples[pair_index][bin_index*2+0];
                obj->frame_complex[bin_index][1] = covs->samples[pair_index][bin_index*2+1];

            }

            fftwf_execute(obj->ifft);

            // Store only the TDoAs in the relevant interval

            memcpy(&obj->cropped_values[obj->tau_max * obj->interpolation_rate], &(obj->frame_real[0]), sizeof(float) * (obj->tau_max * obj->interpolation_rate + 1));
            memcpy(&obj->cropped_values[0], &(obj->frame_real[obj->frame_size * obj->interpolation_rate - obj->tau_max * obj->interpolation_rate]), sizeof(float) * (obj->tau_max * obj->interpolation_rate));

            // Find the TDoA associated to the max cross-correlation amplitude

            max_index = 1;

            for (tau_index = 1; tau_index < (2 * obj->tau_max * obj->interpolation_rate); tau_index++) {
                if (obj->cropped_values[tau_index] > obj->cropped_values[max_index]) {
                    max_index = tau_index;
                }
            }

            // Save the maximum cross-correlation amplitude, and the corresponding TDoA
            // Also save the TDoA before the one with the max cross-correlation amplitude, and the one after
            // (these will be needed later for quadratic interpolation)
            corrs->taus_prev[pair_index] = ((float) (max_index-1) - ((float) obj->tau_max * obj->interpolation_rate)) / ((float) obj->interpolation_rate);
            corrs->ys_prev[pair_index] = obj->cropped_values[max_index-1];
            corrs->taus_max[pair_index] = ((float) (max_index) - ((float) obj->tau_max * obj->interpolation_rate)) / ((float) obj->interpolation_rate);
            corrs->ys_max[pair_index] = obj->cropped_values[max_index];
            corrs->taus_next[pair_index] = ((float) (max_index+1) - ((float) obj->tau_max * obj->interpolation_rate)) / ((float) obj->interpolation_rate);
            corrs->ys_next[pair_index] = obj->cropped_values[max_index+1];

            pair_index++;

        }

    }

    // By default return 0 when returns without error
    return 0;

}

//
// Construct the fcc object
//
// channels_count       Number of channels
// frame_size           Number of samples per frame
// tau_max              Maximum magnitude of the TDoA
//
// (return)             Pointer to the fcc object
//
fcc_obj * fcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max) {

    fcc_obj * obj;
    unsigned int k;
    unsigned int l;

    // Allocate memory for object
    obj = (fcc_obj *) malloc(sizeof(fcc_obj));

    // For now we only support a fixed frame size
    if (frame_size != FCCPHAT_N) {
        printf("Unsupported frame size.\n");
        exit(-1);
    }

    // For now we only support a fixed maximum TDoA
    if (tau_max > FCCPHAT_TAUMAX) {
        printf("Unsupported tau max.\n");
        exit(-1);
    }

    // Copy relevant parameters
    obj->channels_count = channels_count;
    obj->frame_size = frame_size;
    obj->tau_max = FCCPHAT_TAUMAX;

    obj->K = FCCPHAT_K;
    obj->L = FCCPHAT_L;

    // Load bases from constants
    obj->bases = (float **) malloc(sizeof(float *) * obj->K);
    for (k = 0; k < obj->K; k++) {
        obj->bases[k] = (float *) malloc(sizeof(float) * (frame_size/4+1));
        memcpy(obj->bases[k], &(FCCPHAT_BASES[k][0]), sizeof(float) * (frame_size/4+1));
    }

    // Load dictionaries from constants
    obj->dicts = (float **) malloc(sizeof(float *) * obj->L);
    for (l = 0; l < obj->L; l++) {
        obj->dicts[l] = (float *) malloc(sizeof(float) * obj->K * 2);
        memcpy(obj->dicts[l], &(FCCPHAT_DICTS[l][0]), sizeof(float) * obj->K * 2);
    }

    // Return pointer to the object
    return obj;

}

//
// Destroy the fcc object
//
// obj                  Pointer to the fcc object
//
// (return)             None
//
void fcc_destroy(fcc_obj * obj) {

    unsigned int k;
    unsigned int l;

    // Deallocate memory for bases
    for (k = 0; k < obj->K; k++) {
        free((void *) obj->bases[k]);
    }
    free((void *) obj->bases);

    // Deallocate memory for dictionaries
    for (l = 0; l < obj->L; l++) {
        free((void *) obj->dicts[l]);
    }
    free((void *) obj->dicts);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Compute Fast Cross-Correlation
//
// obj                  Pointer to the fcc object
// covs                 Pointer to the covs object that holds the SCM
// corrs                Pointer to the corrs object with the cross-correlation in time-domain
//
// (return)             Returns 0 if no error
//
int fcc_call(fcc_obj * obj, const covs_obj * covs, corrs_obj * corrs) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    float cov_real1;
    float cov_imag1;
    float cov_real2;
    float cov_imag2;

    unsigned int n;
    unsigned int k;
    unsigned int l;

    float x_add_real[(FCCPHAT_N/4+1)];
    float x_add_imag[(FCCPHAT_N/4+1)];
    float x_sub_real[(FCCPHAT_N/4+1)];
    float x_sub_imag[(FCCPHAT_N/4+1)];
    float z_real[FCCPHAT_K];
    float z_imag[FCCPHAT_K];

    float y_real[FCCPHAT_L];
    unsigned int l_max;
    float y_max;

    float current_z_real;
    float current_z_imag;
    float current_y_real;

    pair_index = 0;

    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

        for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {

            // Compute the vectors x_add and x_sub to be used with even and odd bases

            for (n = 0; n < obj->frame_size/4; n++) {

                cov_real1 = covs->samples[pair_index][2*n+0];
                cov_imag1 = covs->samples[pair_index][2*n+1];

                cov_real2 = covs->samples[pair_index][2*(obj->frame_size/2-n)+0];
                cov_imag2 = covs->samples[pair_index][2*(obj->frame_size/2-n)+1];

                x_add_real[n] = cov_real1 + cov_real2;
                x_add_imag[n] = cov_imag1 + cov_imag2;
                x_sub_real[n] = cov_real1 - cov_real2;
                x_sub_imag[n] = cov_imag1 - cov_imag2;

            }

            cov_real1 = covs->samples[pair_index][2*obj->frame_size/4+0];
            cov_imag1 = covs->samples[pair_index][2*obj->frame_size/4+1];

            x_add_real[obj->frame_size/4] = cov_real1;
            x_add_imag[obj->frame_size/4] = cov_imag1;

            x_sub_real[obj->frame_size/4] = cov_real1;
            x_sub_imag[obj->frame_size/4] = cov_imag1;

            // Compute the projection on each even base

            for (k = 0; k < obj->K; k += 2) {
                current_z_real = 0.f;
                current_z_imag = 0.f;

                for (n = 0; n < obj->frame_size/4+1; n++) {
                    current_z_real += x_add_real[n] * obj->bases[k][n];
                    current_z_imag += x_add_imag[n] * obj->bases[k][n];
                }

                z_real[k] = current_z_real;
                z_imag[k] = current_z_imag;
            }

            // Compute the projection on each odd base

            for (k = 1; k < obj->K; k += 2) {
                current_z_real = 0.f;
                current_z_imag = 0.f;

                for (n = 0; n < obj->frame_size/4+1; n++) {
                    current_z_real += -x_sub_imag[n] * obj->bases[k][n];
                    current_z_imag += x_sub_real[n] * obj->bases[k][n];
                }

                z_real[k] = current_z_real;
                z_imag[k] = current_z_imag;
            }

            // Using the vector z, compute the y value for each dictionary element

            for (l = 0; l < obj->L; l++) {
                current_y_real = 0.f;
                for (k = 0; k < obj->K; k++) {
                    current_y_real += z_real[k] * obj->dicts[l][2*k+0] - z_imag[k] * obj->dicts[l][2*k+1];
                }
                y_real[l] = current_y_real;
            }

            // Find the maximum value and corresponding index

            y_max = 0.f;
            l_max = 1;

            for (l = 1; l < obj->L-1; l++) {
                if (y_real[l] > y_max) {
                    y_max = y_real[l];
                    l_max = l;
                }
            }

            // Save the maximum cross-correlation amplitude, and the corresponding TDoA
            // Also save the TDoA before the one with the max cross-correlation amplitude, and the one after
            // (these will be needed later for quadratic interpolation)
            corrs->taus_prev[pair_index] = ((float) (l_max - 1) - (float) (2 * obj->tau_max)) / 2.0;
            corrs->ys_prev[pair_index] = y_real[l_max - 1];
            corrs->taus_max[pair_index] = ((float) l_max - (float) (2 * obj->tau_max)) / 2.0;
            corrs->ys_max[pair_index] = y_real[l_max];
            corrs->taus_next[pair_index] = ((float) (l_max + 1) - (float) (2* obj->tau_max)) / 2.0;
            corrs->ys_next[pair_index] = y_real[l_max + 1];

            pair_index++;

        }

    }

    // By default return 0 when returns without error
    return 0;

}

//
// Construct the quadinterp object
//
// channels_count       Number of channels
//
// (return)             Pointer to the quadinterp object
//
quadinterp_obj * quadinterp_construct(const unsigned int channels_count) {

    quadinterp_obj * obj;

    // Allocate memory for object
    obj = (quadinterp_obj *) malloc(sizeof(quadinterp_obj));

    // Copy relevant parameters
    obj->channels_count = channels_count;

    // Return pointer to object
    return obj;

}

//
// Destroy the quadinterp object
//
// obj                  Pointer to the quadinterp object
//
// (return)             None
//
void quadinterp_destroy(quadinterp_obj * obj) {

    // Deallocate memory for object
    free((void *) obj);

}

//
// Compute the quadratic interpolation
//
// obj                  Pointer to the quadinterp object
// corrs                Pointer to the corrs object with the cross-correlation in time-domain
// taus                 Pointer to the taus object with the tdoas
//
// (return)             Returns 0 if no error
//
int quadinterp_call(quadinterp_obj * obj, const corrs_obj * corrs, taus_obj * taus) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    float y_prev, y_max, y_next;
    float tau_max;
    float delta_tau;
    float tau_hat;

    pair_index = 0;

    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

        for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {

            // Load the maximum value, and neighboor values, and perform quadratic interpolation for each
            // pair of microphones
            y_prev = corrs->ys_prev[pair_index];
            y_max = corrs->ys_max[pair_index];
            y_next = corrs->ys_next[pair_index];
            tau_max = corrs->taus_max[pair_index];
            delta_tau = corrs->taus_next[pair_index] - corrs->taus_max[pair_index];

            tau_hat = tau_max + (delta_tau / 2.0) * (y_prev - y_next) / (y_prev - 2.0 * y_max + y_next);

            taus->taus[pair_index] = tau_hat;
            taus->ys[pair_index] = y_max;

            pair_index++;

        }

    }

    // By default return 0 when returns without error
    return 0;

}

//
// Construct the csv object
//
// file_name            Path of the file name to write the CSV content
// channels_count       Number of channels
//
// (return)             Pointer to the csv object
//
csv_obj * csv_construct(const char * file_name, const unsigned int channels_count) {

    csv_obj * obj;
    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    // Allocate memory for object
    obj = (csv_obj *) malloc(sizeof(csv_obj));

    // Create CSV file
    obj->file_pointer = fopen(file_name, "w");
    if (obj->file_pointer == NULL) {
        printf("Unable to open file %s\n", file_name);
        exit(-1);
    }

    // Copy relevant parameters
    obj->frame_index = 0;
    obj->channels_count = channels_count;

    // Write CSV header to the file
    fprintf(obj->file_pointer, "frame_index, ");

    pair_index = 0;

    for (channel_index1 = 0; channel_index1 < channels_count; channel_index1++) {

        for (channel_index2 = (channel_index1 + 1); channel_index2 < channels_count; channel_index2++) {

            pair_index++;
            fprintf(obj->file_pointer, "tau%03u", pair_index);

            if (pair_index != channels_count * (channels_count-1) / 2) {
                fprintf(obj->file_pointer, ", ");
            }

        }

    }

    fprintf(obj->file_pointer, "\n");

    // Return pointer to object
    return obj;

}

//
// Destroy the csv object
//
// obj                  Pointer to the csv object
//
// (return)             None
//
void csv_destroy(csv_obj * obj) {

    // Close CSV file
    fclose(obj->file_pointer);

    // Deallocate memory for object
    free((void *) obj);

}

//
// Write the TDoAs to the CSV file
//
// obj                  Pointer to the csv object
// taus                 Pointer to the taus object which contains the TDoAs
//
// (return)             Returns 0 if no error
//
int csv_write(csv_obj * obj, taus_obj * taus) {

    unsigned int channel_index1;
    unsigned int channel_index2;
    unsigned int pair_index;

    obj->frame_index++;

    // Write a new line with the frame index first
    fprintf(obj->file_pointer, "%u, ", obj->frame_index);

    // Then write all TDoAs in order on the same line
    pair_index = 0;

    for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

        for (channel_index2 = (channel_index1 + 1); channel_index2 < obj->channels_count; channel_index2++) {

            fprintf(obj->file_pointer, "%+1.4f", taus->taus[pair_index]);

            pair_index++;

            if (pair_index != obj->channels_count * (obj->channels_count-1) / 2) {
                fprintf(obj->file_pointer, ", ");
            }

        }

    }

    fprintf(obj->file_pointer, "\n");

    // By default return 0 when returns without error
    return 0;

}
