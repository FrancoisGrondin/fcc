#include <fccphat/signal.h>

hops_obj * hops_construct(const unsigned int channels_count, const unsigned int hop_size) {

	hops_obj * obj;
	unsigned int channel_index;

	obj = (hops_obj *) malloc(sizeof(hops_obj));

	obj->channels_count = channels_count;
	obj->hop_size = hop_size;

	obj->samples = (float **) malloc(sizeof(float *) * channels_count);
	for (channel_index = 0; channel_index < channels_count; channel_index++) {
		obj->samples[channel_index] = (float *) malloc(sizeof(float) * hop_size);
		memset(obj->samples[channel_index], 0x00, sizeof(float) * hop_size);
	}

	return obj;

}

void hops_destroy(hops_obj * obj) {

	unsigned int channel_index;

	for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
		free((void *) obj->samples[channel_index]);
	}
	free((void *) obj->samples);

	free((void *) obj);

}

void hops_printf(const hops_obj * obj) {

	unsigned int channel_index;
	unsigned int sample_index;

	for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
		
		for (sample_index = 0; sample_index < obj->hop_size; sample_index++){

			printf("[%02u]-(%03u): %1.3f\n", channel_index, sample_index, obj->samples[channel_index][sample_index]);

		}

	}

}


freqs_obj * freqs_construct(const unsigned int channels_count, const unsigned int frame_size) {

	freqs_obj * obj;
	unsigned int channel_index;

	obj = (freqs_obj *) malloc(sizeof(freqs_obj));

	obj->channels_count = channels_count;
	obj->frame_size = frame_size;

	obj->samples = (float **) malloc(sizeof(float *) * channels_count);
	for (channel_index = 0; channel_index < channels_count; channel_index++) {
		obj->samples[channel_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
		memset(obj->samples[channel_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
	}

	return obj;

}

void freqs_destroy(freqs_obj * obj) {

	unsigned int channel_index;

	for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
		free((void *) obj->samples[channel_index]);
	}
	free((void *) obj->samples);

	free((void *) obj);

}

void freqs_printf(const freqs_obj * obj) {

	unsigned int channel_index;
	unsigned int bin_index;

	for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
		
		for (bin_index = 0; bin_index < obj->frame_size/2+1; bin_index++){

			printf("[%02u]-(%03u): (%+1.3f , %+1.3f)\n", channel_index, bin_index, obj->samples[channel_index][bin_index*2+0], obj->samples[channel_index][bin_index*2+1]);

		}

	}

}


covs_obj * covs_construct(const unsigned int channels_count, const unsigned int frame_size) {

	covs_obj * obj;
	unsigned int channel_index1;
	unsigned int channel_index2;
	unsigned int pair_index;

	obj = (covs_obj *) malloc(sizeof(covs_obj));

	obj->channels_count = channels_count;
	obj->frame_size = frame_size;

	pair_index = 0;
	obj->samples = (float **) malloc(sizeof(float *) * channels_count * (channels_count-1) / 2);
	for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
		for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
			obj->samples[pair_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
			memset(obj->samples[pair_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
			pair_index++;
		}
	}

	return obj;

}

void covs_destroy(covs_obj * obj) {

	unsigned int channel_index1;
	unsigned int channel_index2;
	unsigned int pair_index;

	pair_index = 0;
	for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {
		for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {
			free((void *) obj->samples[pair_index]);
			pair_index++;
		}
	}
	free((void *) obj->samples);

	free((void *) obj);

}

void covs_printf(const covs_obj * obj) {

	unsigned int pair_index;
	unsigned int bin_index;

	for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count-1)/2); pair_index++) {
		
		for (bin_index = 0; bin_index < obj->frame_size/2+1; bin_index++){

			printf("[%02u]-(%03u): (%+1.3f , %+1.3f)\n", pair_index, bin_index, obj->samples[pair_index][bin_index*2+0], obj->samples[pair_index][bin_index*2+1]);

		}

	}

}

corrs_obj * corrs_construct(const unsigned int channels_count) {

	corrs_obj * obj;

	obj = (corrs_obj *) malloc(sizeof(corrs_obj));

	obj->channels_count = channels_count;

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

	return obj;

}

void corrs_destroy(corrs_obj * obj) {

	free((void *) obj->taus_prev);
	free((void *) obj->ys_prev);
	free((void *) obj->taus_max);
	free((void *) obj->ys_max);
	free((void *) obj->taus_next);
	free((void *) obj->ys_next);

	free((void *) obj);

}

void corrs_printf(const corrs_obj * obj) {

	unsigned int pair_index;

	for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count-1)/2); pair_index++) {

		printf("[%02u] taus_prev=%+1.3f ys_prev=%+1.3f taus_max=%+1.3f ys_max=%+1.3f taus_next=%+1.3f ys_next=%+1.3f\n", 
			pair_index, obj->taus_prev[pair_index], obj->ys_prev[pair_index], obj->taus_max[pair_index], obj->ys_max[pair_index], obj->taus_next[pair_index], obj->ys_next[pair_index]);

	}

}


taus_obj * taus_construct(const unsigned int channels_count) {

	taus_obj * obj;

	obj = (taus_obj *) malloc(sizeof(taus_obj));

	obj->channels_count = channels_count;

	obj->taus = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
	memset(obj->taus, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);
	obj->ys = (float *) malloc(sizeof(float) * channels_count * (channels_count-1) / 2);
	memset(obj->ys, 0x00, sizeof(float) * channels_count * (channels_count-1) / 2);

	return obj;

}

void taus_destroy(taus_obj * obj) {

	free((void *) obj->taus);
	free((void *) obj->ys);

	free((void *) obj);

}

void taus_printf(const taus_obj * obj) {

	unsigned int pair_index;

	for (pair_index = 0; pair_index < (obj->channels_count * (obj->channels_count - 1)/2); pair_index++) {

		printf("[%02u] taus=%1.3f ys=%1.3f\n", pair_index, obj->taus[pair_index], obj->ys[pair_index]);

	}

}