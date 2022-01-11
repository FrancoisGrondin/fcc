#ifndef __FCCPHAT_SYSTEM
#define __FCCPHAT_SYSTEM

	#include <stdlib.h>
	#include <stdio.h>
	#include <string.h>
	#include <math.h>
	#include <fftw3.h>
	#include "signal.h"
	#include "const.h"

	typedef struct wav_header {

		char chunk_id[4];
		unsigned int chunk_size;
		char format[4];
		char subchunk1_id[4];
		unsigned int subchunk1_size;
		unsigned short audio_format;
		unsigned short num_channels;
		unsigned int sample_rate;
		unsigned int byte_rate;
		unsigned short block_align;
		unsigned short bits_per_sample;
		char subchunk2_id[4];
		unsigned int subchunk2_size;

	} wav_header;

	typedef struct wav_obj {

		unsigned int sample_rate;
		unsigned int num_channels;
		unsigned int bits_per_sample;
		unsigned int hop_size;

		FILE * file_pointer;
		short * buffer;

	} wav_obj;
	
	typedef struct stft_obj {

		unsigned int channels_count;
		unsigned int frame_size;
		unsigned int hop_size;

		float * window;
		float ** frames;

        fftwf_plan fft;
        float * frame_real;
        fftwf_complex * frame_complex;

	} stft_obj;	

	typedef struct scm_obj {

		unsigned int channels_count;
		unsigned int frame_size;
		float alpha;

		float ** cross_spectrum;

	} scm_obj;

	typedef struct phat_obj {

		unsigned int channels_count;
		unsigned int frame_size;

	} phat_obj;

	typedef struct gcc_obj {

		unsigned int channels_count;
		unsigned int frame_size;
		unsigned int tau_max;
		unsigned int interpolation_rate;

        fftwf_plan ifft;
        fftwf_complex * frame_complex;
        float * frame_real;

        float * cropped_values;

	} gcc_obj;

	typedef struct fcc_obj {

		unsigned int channels_count;
		unsigned int frame_size;
		unsigned int tau_max;

		unsigned int K;
		unsigned int L;
		float ** bases;
		float ** dicts;

		float * x_add;
		float * x_sub;
		float * z;
		float * y;

	} fcc_obj;

	typedef struct quadinterp_obj {

		unsigned int channels_count;

	} quadinterp_obj;

	wav_obj * wav_construct(const char * file_name, const unsigned int hop_size);

	void wav_destroy(wav_obj * obj);

	int wav_read(wav_obj * obj, hops_obj * hops);

	stft_obj * stft_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int hop_size);

	void stft_destroy(stft_obj * obj);

	int stft_call(stft_obj * obj, const hops_obj * hops, freqs_obj * freqs);

	scm_obj * scm_construct(const unsigned int channels_count, const unsigned int frame_size, const float alpha);

	void scm_destroy(scm_obj * obj);

	int scm_call(scm_obj * obj, const freqs_obj * freqs, covs_obj * covs);

	phat_obj * phat_construct(const unsigned int channels_count, const unsigned int frame_size);

	void phat_destroy(phat_obj * obj);

	int phat_call(phat_obj * obj, covs_obj * covs, covs_obj * covs_normalized);

	gcc_obj * gcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max, const unsigned int interpolation_rate);

	void gcc_destroy(gcc_obj * obj);

	int gcc_call(gcc_obj * obj, const covs_obj * covs, corrs_obj * corrs);

	fcc_obj * fcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max);

	void fcc_destroy(fcc_obj * obj);

	int fcc_call(fcc_obj * obj, const covs_obj * covs, corrs_obj * corrs);

	quadinterp_obj * quadinterp_construct(const unsigned int channels_count);

	void quadinterp_destroy(quadinterp_obj * obj);

	int quadinterp_call(quadinterp_obj * obj, const corrs_obj * corrs, taus_obj * taus);

#endif