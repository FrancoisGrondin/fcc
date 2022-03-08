#include <fccphat/signal.h>
#include <fccphat/system.h>

#include <unistd.h>
#include <string.h>

int main(int argc, char * argv[]) {

    int opt;
    char step;
    unsigned int repeat;

    stft_obj * stft;
    scmphat_obj * scmphat;
    gcc_obj * gcc;
    fcc_obj * fcc;
    quadinterp_obj * quadinterp;

    hops_obj * hops;
    freqs_obj * freqs;
    covs_obj * covs;
    corrs_obj * corrs;
    taus_obj * taus;

    const unsigned int hop_size = 256;
    const unsigned int frame_size = 512;
    const unsigned int tau_max = 8;
    const unsigned int interpolation_rate = 2;
    const float alpha = 0.1;

    unsigned int r;
    unsigned int channels_count;

    step = 0;
    repeat = 0;
    channels_count = 2;

    while((opt = getopt(argc, argv, "b:c:r:")) != -1) 
    { 
        switch(opt) 
        { 

            case 'b':

                if (strcmp(optarg, "stft") == 0) { step = 1; }
                if (strcmp(optarg, "scmphat") == 0) { step = 2; }
                if (strcmp(optarg, "gcc") == 0) { step = 3; }
                if (strcmp(optarg, "fcc") == 0) { step = 4; }
                if (strcmp(optarg, "qi") == 0) { step = 5; }
                break;

            case 'c':

                channels_count = atoi(optarg);

            case 'r':

                repeat = atoi(optarg);
                break;

        } 
    } 

    stft = stft_construct(channels_count, frame_size, hop_size);
    scmphat = scmphat_construct(channels_count, frame_size, alpha);
    gcc = gcc_construct(channels_count, frame_size, tau_max, interpolation_rate);
    fcc = fcc_construct(channels_count, frame_size, tau_max);
    quadinterp = quadinterp_construct(channels_count);

    hops = hops_construct(channels_count, hop_size);
    freqs = freqs_construct(channels_count, frame_size);
    covs = covs_construct(channels_count, frame_size);
    corrs = corrs_construct(channels_count);
    taus = taus_construct(channels_count);

    for (r = 0; r < repeat; r++) {

        if (step == 1) { stft_call(stft, hops, freqs); }
        if (step == 2) { scmphat_call(scmphat, freqs, covs); }
        if (step == 3) { gcc_call(gcc, covs, corrs); }
        if (step == 4) { fcc_call(fcc, covs, corrs); }
        if (step == 5) { quadinterp_call(quadinterp, corrs, taus); }

    }

    stft_destroy(stft);
    scmphat_destroy(scmphat);
    gcc_destroy(gcc);
    fcc_destroy(fcc);
    quadinterp_destroy(quadinterp);

    hops_destroy(hops);
    freqs_destroy(freqs);
    covs_destroy(covs);
    corrs_destroy(corrs);
    taus_destroy(taus);

    return 0;

}