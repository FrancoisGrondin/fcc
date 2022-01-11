#include <fccphat/signal.h>
#include <fccphat/system.h>

#include <unistd.h>
#include <string.h>

int main(int argc, char * argv[]) {

    int opt;
    char method;
    char filename[1024];
    unsigned int repeat;

    wav_obj * wav;
    stft_obj * stft;
    scm_obj * scm;
    phat_obj * phat;
    gcc_obj * gcc;
    fcc_obj * fcc;
    quadinterp_obj * quadinterp;

    hops_obj * hops;
    freqs_obj * freqs;
    covs_obj * covs;
    covs_obj * covs_normalized;
    corrs_obj * corrs;
    taus_obj * taus;

    const unsigned int hop_size = 256;
    const unsigned int frame_size = 512;
    const unsigned int channels_count = 2;
    const unsigned int sample_rate = 16000;
    const unsigned int tau_max = 8;
    const unsigned int interpolation_rate = 2;
    const float alpha = 1.0;

    unsigned int r;

    method = 0x00;
    repeat = 1;

    while((opt = getopt(argc, argv, "i:m:r:")) != -1) 
    { 
        switch(opt) 
        { 
            
            case 'i':

                strcpy(filename, optarg);
                break;

            case 'm':

                if (strcmp(optarg, "gcc") == 0) { method = 'g'; }
                if (strcmp(optarg, "fcc") == 0) { method = 'f'; }
                break;

            case 'r':

                repeat = atoi(optarg);
                break;

        } 
    } 

    if (method == 0x00) {
        printf("Invalid method. Choose either \"fcc\" or \"gcc\"\n");
        exit(-1);        
    }

    if( access( filename, F_OK ) != 0 ) {
        printf("File %s does not exists\n", filename);
        exit(-1);
    }

    wav = wav_construct(filename, hop_size);

    if (wav->sample_rate != sample_rate) {
        printf("Invalid sample rate (was expecting %u, got %u).\n", sample_rate, wav->sample_rate);
    }
    if (wav->num_channels != channels_count) {
        printf("Invalid number of channels (was expecting %u, got %u).\n", channels_count, wav->num_channels);
    }

    stft = stft_construct(channels_count, frame_size, hop_size);
    scm = scm_construct(channels_count, frame_size, alpha);
    phat = phat_construct(channels_count, frame_size);
    gcc = gcc_construct(channels_count, frame_size, tau_max, interpolation_rate);
    fcc = fcc_construct(channels_count, frame_size, tau_max);
    quadinterp = quadinterp_construct(channels_count);

    hops = hops_construct(channels_count, hop_size);
    freqs = freqs_construct(channels_count, frame_size);
    covs = covs_construct(channels_count, frame_size);
    covs_normalized = covs_construct(channels_count, frame_size);
    corrs = corrs_construct(channels_count);
    taus = taus_construct(channels_count);

    while (wav_read(wav, hops) == 0) {

        stft_call(stft, hops, freqs);
        scm_call(scm, freqs, covs);
        phat_call(phat, covs, covs_normalized);

        for (r = 0; r < repeat; r++) {

            if (method == 'g') {
                gcc_call(gcc, covs_normalized, corrs);
            }
            if (method == 'f') {
                fcc_call(fcc, covs_normalized, corrs);
            }
            
        }

        quadinterp_call(quadinterp, corrs, taus);

    }

    stft_destroy(stft);
    scm_destroy(scm);
    phat_destroy(phat);
    gcc_destroy(gcc);
    fcc_destroy(fcc);
    quadinterp_destroy(quadinterp);

    hops_destroy(hops);
    freqs_destroy(freqs);
    covs_destroy(covs);
    covs_destroy(covs_normalized);
    corrs_destroy(corrs);
    taus_destroy(taus);

    return 0;

}