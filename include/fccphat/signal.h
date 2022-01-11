#ifndef __FCCPHAT_SIGNAL
#define __FCCPHAT_SIGNAL

	#include <stdlib.h>
	#include <string.h>

	typedef struct hops_obj {

		unsigned int channels_count;
		unsigned int hop_size;
		float ** samples;

	} hops_obj;

	typedef struct freqs_obj {

		unsigned int channels_count;
		unsigned int frame_size;
		float ** samples;

	} freqs_obj;

	typedef struct covs_obj {

		unsigned int channels_count;
		unsigned int frame_size;
		float ** samples;

	} covs_obj;

	typedef struct corrs_obj {

		unsigned int channels_count;

		float * taus_prev;
		float * ys_prev;

		float * taus_max;
		float * ys_max;

		float * taus_next;
		float * ys_next;

	} corrs_obj;

	typedef struct taus_obj {

		unsigned int channels_count;

		float * taus;
		float * ys;

	} taus_obj;

	hops_obj * hops_construct(const unsigned int channels_count, const unsigned int hop_size);

	void hops_destroy(hops_obj * obj);

	freqs_obj * freqs_construct(const unsigned int channels_count, const unsigned int frame_size);

	void freqs_destroy(freqs_obj * obj);

	covs_obj * covs_construct(const unsigned int channels_count, const unsigned int frame_size);

	void covs_destroy(covs_obj * obj);

	corrs_obj * corrs_construct(const unsigned int channels_count);

	void corrs_destroy(corrs_obj * obj);

	taus_obj * taus_construct(const unsigned int channels_count);

	void taus_destroy(taus_obj * obj);

#endif