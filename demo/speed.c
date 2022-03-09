#include <fccphat/signal.h>
#include <fccphat/system.h>

#include <unistd.h>
#include <string.h>
#include <time.h>

int main(int argc, char * argv[]) {

    int opt;
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
    const float alpha = 0.1f;

    unsigned int r;
    unsigned int channels_count;

    clock_t start;
    float stft_duration;
    float scmphat_duration;
    float gcc_duration;
    float fcc_duration;
    float quadinterp_duration;

    repeat = 0;
    channels_count = 2;

    while((opt = getopt(argc, argv, "c:r:")) != -1)
    {
        switch(opt)
        {

            case 'c':

                channels_count = atoi(optarg);
                break;

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

    start = clock();
    for (r = 0; r < repeat; r++) {
        stft_call(stft, hops, freqs);
    }
    stft_duration = (float)(clock() - start) / CLOCKS_PER_SEC;

    start = clock();
    for (r = 0; r < repeat; r++) {
        scmphat_call(scmphat, freqs, covs);
    }
    scmphat_duration = (float)(clock() - start) / CLOCKS_PER_SEC;

    start = clock();
    for (r = 0; r < repeat; r++) {
        gcc_call(gcc, covs, corrs);
    }
    gcc_duration = (float)(clock() - start) / CLOCKS_PER_SEC;

    start = clock();
    for (r = 0; r < repeat; r++) {
        fcc_call(fcc, covs, corrs);
    }
    fcc_duration = (float)(clock() - start) / CLOCKS_PER_SEC;

    start = clock();
    for (r = 0; r < repeat; r++) {
        quadinterp_call(quadinterp, corrs, taus);
    }
    quadinterp_duration = (float)(clock() - start) / CLOCKS_PER_SEC;

    printf("***** Duration *****\n");
    printf("    stft:       %1.3f s\n", stft_duration);
    printf("    scmphat:    %1.3f s\n", scmphat_duration);
    printf("    gcc:        %1.3f s\n", gcc_duration);
    printf("    fcc:        %1.3f s\n", fcc_duration);
    printf("    quadinterp: %1.3f s\n\n", quadinterp_duration);

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
