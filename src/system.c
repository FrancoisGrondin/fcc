#include <fccphat/system.h>

wav_obj * wav_construct(const char * file_name, const unsigned int hop_size) {

	wav_obj * obj;
	wav_header hdr;
	size_t rtn;

	obj = (wav_obj *) malloc(sizeof(wav_obj));

	obj->file_pointer = fopen(file_name, "rb");

	rtn = fread(&hdr, sizeof(char), sizeof(wav_header), obj->file_pointer);

	obj->sample_rate = hdr.sample_rate;
	obj->num_channels = hdr.num_channels;
	obj->bits_per_sample = hdr.bits_per_sample;
	obj->hop_size = hop_size;

	if (obj->bits_per_sample != 16) {
		printf("Unsupported PCM format.\n");
		exit(-1);
	}

	obj->buffer = (short *) malloc(sizeof(short) * obj->hop_size * obj->num_channels);

	return obj;

}

void wav_destroy(wav_obj * obj) {

	fclose(obj->file_pointer);
	free((void *) obj->buffer);

	free((void *) obj);

}

int wav_read(wav_obj * obj, hops_obj * hops) {
	
	size_t rtn;
	unsigned int channel_index;
	unsigned int sample_index;

	rtn = fread(obj->buffer, obj->bits_per_sample/8, obj->hop_size * obj->num_channels, obj->file_pointer);

	for (sample_index = 0; sample_index < obj->hop_size; sample_index++) {
		for (channel_index = 0; channel_index < obj->num_channels; channel_index++) {
			hops->samples[channel_index][sample_index] = ((float) obj->buffer[sample_index * obj->num_channels + channel_index]) / 32768.0;
		}
	}

	return (rtn == obj->hop_size * obj->num_channels) ? 0:-1;

}

stft_obj * stft_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int hop_size) {

	stft_obj * obj;
	unsigned int channel_index;
	unsigned int sample_index;

	obj = (stft_obj *) malloc(sizeof(stft_obj));

	obj->channels_count = channels_count;
	obj->frame_size = frame_size;
	obj->hop_size = hop_size;

	obj->window = (float *) malloc(sizeof(float) * frame_size);
	for (sample_index = 0; sample_index < frame_size; sample_index++) {
		obj->window[sample_index] = powf(sinf((float) sample_index * M_PI / ((float) (frame_size - 1))), 2.0);
	}

	obj->frames = (float **) malloc(sizeof(float *) * channels_count);
	for (channel_index = 0; channel_index < channels_count; channel_index++) {
		obj->frames[channel_index] = (float *) malloc(sizeof(float) * frame_size);
		memset((void *) obj->frames[channel_index], 0x00, sizeof(float) * frame_size);
	}

	obj->frame_real = (float *) fftwf_malloc(sizeof(float) * frame_size);
	obj->frame_complex = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * (frame_size/2+1));
	obj->fft = fftwf_plan_dft_r2c_1d(frame_size, obj->frame_real, obj->frame_complex, FFTW_ESTIMATE);

	return obj;

}

void stft_destroy(stft_obj * obj) {

	unsigned int channel_index;

	free((void *) obj->window);

	for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
		free((void *) obj->frames[channel_index]);
	}
	free((void *) obj->frames);

    fftwf_free(obj->frame_real);
    fftwf_free(obj->frame_complex);
    fftwf_destroy_plan(obj->fft);
    fftwf_cleanup();

	free((void *) obj);

}

int stft_call(stft_obj * obj, const hops_obj * hops, freqs_obj * freqs) {

	unsigned int channel_index;
	unsigned int sample_index;
	unsigned int bin_index;

	for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {

		memmove(&(obj->frames[channel_index][0]), &(obj->frames[channel_index][obj->hop_size]), sizeof(float) * (obj->frame_size - obj->hop_size));
		memcpy(&(obj->frames[channel_index][obj->frame_size - obj->hop_size]), &(hops->samples[channel_index][0]), sizeof(float) * obj->hop_size);

		for (sample_index = 0; sample_index < obj->frame_size; sample_index++) {
			obj->frame_real[sample_index] = obj->window[sample_index] * obj->frames[channel_index][sample_index];
		}

		fftwf_execute(obj->fft);

		for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {
			freqs->samples[channel_index][bin_index * 2 + 0] = obj->frame_complex[bin_index][0];
			freqs->samples[channel_index][bin_index * 2 + 1] = obj->frame_complex[bin_index][1];
		}

	}

	return 0;

}

scmphat_obj * scmphat_construct(const unsigned int channels_count, const unsigned int frame_size, const float alpha) {

	scmphat_obj * obj;
	unsigned int channel_index1;
	unsigned int channel_index2;
	unsigned int pair_index;
	unsigned int channel_index;

	obj = (scmphat_obj *) malloc(sizeof(scmphat_obj));

	obj->channels_count = channels_count;
	obj->frame_size = frame_size;
	obj->alpha = alpha;

	obj->method = 'p';
	if (alpha == 1.0f) {
		if (channels_count > 2) {
			obj->method = 'c';
		}
	}

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
	if (obj->method == 'c') {

		obj->cross_spectrum = (float **) malloc(sizeof(float *) * channels_count);
		for (channel_index = 0; channel_index < channels_count; channel_index++) {
			obj->cross_spectrum[channel_index] = (float *) malloc(sizeof(float) * (frame_size/2+1) * 2);
			memset((void *) obj->cross_spectrum[channel_index], 0x00, sizeof(float) * (frame_size/2+1) * 2);
		}

	}

	return obj;

}

void scmphat_destroy(scmphat_obj * obj) {

	unsigned int channel_index1;
	unsigned int channel_index2;
	unsigned int pair_index;
	unsigned int channel_index;

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
	if (obj->method == 'c') {

		for (channel_index = 0; channel_index < obj->channels_count; channel_index++) {
			free((void *) obj->cross_spectrum[channel_index]);
		}
		free((void *) obj->cross_spectrum);

	}

	free((void *) obj);

}

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

	if (obj->method == 'p') {

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
	if (obj->method == 'c') {

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

	return 0;

}

gcc_obj * gcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max, const unsigned int interpolation_rate) {

	gcc_obj * obj;

	obj = (gcc_obj *) malloc(sizeof(gcc_obj));

	obj->channels_count = channels_count;
	obj->frame_size = frame_size;
	obj->tau_max = tau_max;
	obj->interpolation_rate = interpolation_rate;

	obj->frame_complex = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * (frame_size/2+1) * interpolation_rate);
	obj->frame_real = (float *) fftwf_malloc(sizeof(float) * frame_size * interpolation_rate);
	obj->ifft = fftwf_plan_dft_c2r_1d(frame_size * interpolation_rate, obj->frame_complex, obj->frame_real, FFTW_ESTIMATE);	

	obj->cropped_values = (float *) malloc(sizeof(float) * (tau_max * 2 * interpolation_rate + 1));

	return obj;

}

void gcc_destroy(gcc_obj * obj) {

	free((void *) obj->frame_complex);
	free((void *) obj->frame_real);
    fftwf_destroy_plan(obj->ifft);
    fftwf_cleanup();

    free((void *) obj->cropped_values);

	free((void *) obj);

}

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

			memset(obj->frame_complex, 0x00, sizeof(fftwf_complex) * (obj->frame_size/2+1) * obj->interpolation_rate);

			for (bin_index = 0; bin_index < (obj->frame_size/2+1); bin_index++) {

				obj->frame_complex[bin_index][0] = covs->samples[pair_index][bin_index*2+0];
				obj->frame_complex[bin_index][1] = covs->samples[pair_index][bin_index*2+1];		

			}

			fftwf_execute(obj->ifft);

			memcpy(&obj->cropped_values[obj->tau_max * obj->interpolation_rate], &(obj->frame_real[0]), sizeof(float) * (obj->tau_max * obj->interpolation_rate + 1));
			memcpy(&obj->cropped_values[0], &(obj->frame_real[obj->frame_size * obj->interpolation_rate - obj->tau_max * obj->interpolation_rate]), sizeof(float) * (obj->tau_max * obj->interpolation_rate));

			max_index = 1;

			for (tau_index = 1; tau_index < (2 * obj->tau_max * obj->interpolation_rate); tau_index++) {
				if (obj->cropped_values[tau_index] > obj->cropped_values[max_index]) {
					max_index = tau_index;
				}
			}

			corrs->taus_prev[pair_index] = ((float) (max_index-1) - ((float) obj->tau_max * obj->interpolation_rate)) / ((float) obj->interpolation_rate);
			corrs->ys_prev[pair_index] = obj->cropped_values[max_index-1];
			corrs->taus_max[pair_index] = ((float) (max_index) - ((float) obj->tau_max * obj->interpolation_rate)) / ((float) obj->interpolation_rate);
			corrs->ys_max[pair_index] = obj->cropped_values[max_index];
			corrs->taus_next[pair_index] = ((float) (max_index+1) - ((float) obj->tau_max * obj->interpolation_rate)) / ((float) obj->interpolation_rate);
			corrs->ys_next[pair_index] = obj->cropped_values[max_index+1];

			pair_index++;

		}

	}

	return 0;

}

fcc_obj * fcc_construct(const unsigned int channels_count, const unsigned int frame_size, const unsigned int tau_max) {

	fcc_obj * obj;
	unsigned int k;
	unsigned int l;

	obj = (fcc_obj *) malloc(sizeof(fcc_obj));

	if (frame_size != FCCPHAT_N) {
		printf("Unsupported frame size.\n");
		exit(-1);
	}

	if (tau_max > FCCPHAT_TAUMAX) {
		printf("Unsupported tau max.\n");
		exit(-1);
	}

	obj->channels_count = channels_count;
	obj->frame_size = frame_size;
	obj->tau_max = FCCPHAT_TAUMAX;

	obj->K = FCCPHAT_K;
	obj->L = FCCPHAT_L;

	obj->bases = (float **) malloc(sizeof(float *) * obj->K);
	for (k = 0; k < obj->K; k++) {
		obj->bases[k] = (float *) malloc(sizeof(float) * (frame_size/4+1));
		memcpy(obj->bases[k], &(FCCPHAT_BASES[k][0]), sizeof(float) * (frame_size/4+1));
	}

	obj->dicts = (float **) malloc(sizeof(float *) * obj->L);
	for (l = 0; l < obj->L; l++) {
		obj->dicts[l] = (float *) malloc(sizeof(float) * obj->K * 2);
		memcpy(obj->dicts[l], &(FCCPHAT_DICTS[l][0]), sizeof(float) * obj->K * 2);
	}

	return obj;

}

void fcc_destroy(fcc_obj * obj) {

	unsigned int k;
	unsigned int l;

	for (k = 0; k < obj->K; k++) {
		free((void *) obj->bases[k]);
	}
	free((void *) obj->bases);

	for (l = 0; l < obj->L; l++) {
		free((void *) obj->dicts[l]);
	}
	free((void *) obj->dicts);

	free((void *) obj);

}

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

	pair_index = 0;

	for (channel_index1 = 0; channel_index1 < obj->channels_count; channel_index1++) {

		for (channel_index2 = (channel_index1+1); channel_index2 < obj->channels_count; channel_index2++) {

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

			cov_real1 = covs->samples[pair_index][2*obj->frame_size/2+0];
			cov_imag1 = covs->samples[pair_index][2*obj->frame_size/2+1];

			x_add_real[obj->frame_size/4] = cov_real1;
			x_add_imag[obj->frame_size/4] = cov_imag1;

			x_sub_real[obj->frame_size/4] = cov_real1;
			x_sub_imag[obj->frame_size/4] = cov_imag1;

			for (k = 0; k < obj->K; k++) {

				z_real[k] = 0.0;
				z_imag[k] = 0.0;

				if (k % 2 == 0) {
					for (n = 0; n < obj->frame_size/4+1; n++) {
						z_real[k] += x_add_real[n] * obj->bases[k][n];
						z_imag[k] += x_add_imag[n] * obj->bases[k][n];
					}
				}
				else {
					for (n = 0; n < obj->frame_size/4+1; n++) {
						z_real[k] += x_sub_real[n] * obj->bases[k][n];
						z_imag[k] += x_sub_imag[n] * obj->bases[k][n];
					}					
				}	


			}

			for (l = 0; l < obj->L; l++) {
				y_real[l] = 0.0;
				for (k = 0; k < obj->K; k++) {
					y_real[l] += z_real[k] * obj->dicts[l][2*k+0] - z_imag[k] * obj->dicts[l][2*k+1];
				}
			}

			y_max = 0.0;
			l_max = 1;

			for (l = 1; l < obj->L-1; l++) {
				if (y_real[l] > y_max) {
					y_max = y_real[l];
					l_max = l;
				}
			}	

			corrs->taus_prev[pair_index] = ((float) (l_max - 1) - (float) (2 * obj->tau_max)) / 2.0;
			corrs->ys_prev[pair_index] = y_real[l_max - 1];
			corrs->taus_max[pair_index] = ((float) l_max - (float) (2 * obj->tau_max)) / 2.0;
			corrs->ys_max[pair_index] = y_real[l_max];
			corrs->taus_next[pair_index] = ((float) (l_max + 1) - (float) (2* obj->tau_max)) / 2.0;
			corrs->ys_next[pair_index] = y_real[l_max + 1];

			pair_index++;

		}

	}

	return 0;

}

quadinterp_obj * quadinterp_construct(const unsigned int channels_count) {

	quadinterp_obj * obj;

	obj = (quadinterp_obj *) malloc(sizeof(quadinterp_obj));

	obj->channels_count = channels_count;

	return obj;

}

void quadinterp_destroy(quadinterp_obj * obj) {

	free((void *) obj);

}

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

			y_prev = corrs->ys_prev[pair_index];
			y_max = corrs->ys_max[pair_index];
			y_next = corrs->ys_next[pair_index];
			tau_max = corrs->taus_max[pair_index];
			delta_tau = corrs->taus_next[pair_index] - corrs->taus_max[pair_index];

			tau_hat = tau_max + (delta_tau / 2.0) * (y_prev - y_next) / (y_prev - 2.0 * y_max + y_next);

			taus->taus[pair_index] = tau_hat;

			pair_index++;

		}

	}

	return 0;

}

csv_obj * csv_construct(const char * file_name, const unsigned int channels_count) {

	csv_obj * obj;
	unsigned int channel_index1;
	unsigned int channel_index2;
	unsigned int pair_index;

	obj = (csv_obj *) malloc(sizeof(csv_obj));

	obj->file_pointer = fopen(file_name, "w");
	obj->frame_index = 0;
	obj->channels_count = channels_count;

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

	return obj;

}

void csv_destroy(csv_obj * obj) {

	fclose(obj->file_pointer);

	free((void *) obj);

}

int csv_write(csv_obj * obj, taus_obj * taus) {

	unsigned int channel_index1;
	unsigned int channel_index2;
	unsigned int pair_index;

	obj->frame_index++;

	fprintf(obj->file_pointer, "%u, ", obj->frame_index);

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

	return 0;

}